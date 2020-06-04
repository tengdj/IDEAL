/*
 * Parser.cpp
 *
 *  Created on: May 9, 2020
 *      Author: teng
 */

#include "../geometry/MyPolygon.h"
#include <fstream>
#include "../index/RTree.h"
#include <queue>
#include <boost/program_options.hpp>
#include "cuda/mygpu.h"

namespace po = boost::program_options;
using namespace std;

#define MAX_THREAD_NUM 100

int element_size = 100;

// some shared parameters
pthread_mutex_t poly_lock;
bool stop = false;

pthread_mutex_t report_lock;
long query_count = 0;

vector<gpu_info *> gpus;
pthread_mutex_t gpu_lock;

class queue_element{
public:
	int id = 0;
	vector<MyPolygon *> polys;
	queue_element(int i){
		id = i;
	}
	~queue_element(){
		for(MyPolygon *p:polys){
			delete p;
		}
		polys.clear();
	}
};
queue<queue_element *> shared_queue;

RTree<MyPolygon *, double, 2, double> tree;
void contain_batch_gpu(gpu_info *gpu, double *data, uint *offset_size, int *result, size_t total_vertice_num, int pair_num);
void load_source_togpu(gpu_info *gpu, vector<MyPolygon *> &source);

int process_with_gpu(query_context *ctx){

	int pair_num = ctx->candidates.size();
	if(pair_num==0){
		return 0;
	}
	int *result = new int[pair_num];
	uint *offset_size = new uint[pair_num*4];
	uint dataoffset = 0;
	uint total_num_vertices = 0;
	for(int i=0;i<pair_num;i++){
		if(i==0||ctx->candidates[i].second->getid()!=ctx->candidates[i-1].second->getid()){
			total_num_vertices += ctx->candidates[i].second->boundary->num_vertices;
		}
	}

	double *tmpdata = new double[total_num_vertices*2];
	for(int i=0;i<pair_num;i++){
		offset_size[i*4] = ctx->candidates[i].first->offset;
		offset_size[i*4+1] = ctx->candidates[i].first->boundary->num_vertices;
		offset_size[i*4+3] = ctx->candidates[i].second->boundary->num_vertices;
		if(i==0||ctx->candidates[i].second->getid()!=ctx->candidates[i-1].second->getid()){
			offset_size[i*4+2] = dataoffset;
			int num_vertices = ctx->candidates[i].second->boundary->num_vertices;
			memcpy((char *)(tmpdata+dataoffset), (char *)(ctx->candidates[i].second->boundary->x), num_vertices*sizeof(double));
			dataoffset += num_vertices;
			memcpy((char *)(tmpdata+dataoffset), (char *)(ctx->candidates[i].second->boundary->y), num_vertices*sizeof(double));
			dataoffset += num_vertices;
		}else{
			offset_size[i*4+2] = dataoffset-offset_size[i*4+3]*2;
		}
	}
	assert(dataoffset==total_num_vertices*2);
	ctx->candidates.clear();

	int found = 0;
	pthread_mutex_lock(&gpu_lock);
	contain_batch_gpu(gpus[0],tmpdata,offset_size,result,total_num_vertices,pair_num);
	pthread_mutex_unlock(&gpu_lock);
	for(int i=0;i<pair_num;i++){
		found += result[i];
	}
	delete []result;
	delete []offset_size;
	delete []tmpdata;
	return found;
}

bool MySearchCallback(MyPolygon *poly, void* arg){
	query_context *ctx = (query_context *)arg;
	// MBB check
	if(!poly->getMBB()->contain(ctx->target->getMBB())){
		return true;
	}
	// query with parition
	if(ctx->use_partition){
		if(!poly->ispartitioned()){
			poly->partition(ctx->vpr);
		}
		bool find = poly->contain_try_partition(ctx->target, ctx);
		if(ctx->partition_determined){
			ctx->found += find;
			return true;
		}
	}

	if(!ctx->gpu){
		ctx->found += poly->contain(ctx->target, ctx);
	}else{
		ctx->candidates.push_back(std::make_pair(poly, ctx->target));
	}
	// keep going until all hit objects are found
	return true;
}

void *query(void *args){
	query_context *ctx = (query_context *)args;
	pthread_mutex_lock(&poly_lock);
	log("thread %d is started",ctx->thread_id);
	pthread_mutex_unlock(&poly_lock);
	vector<queue_element *> elems;
	int local_count = 0;
	while(!stop||shared_queue.size()>0){
		queue_element *elem = NULL;
		pthread_mutex_lock(&poly_lock);
		if(!shared_queue.empty()){
			elem = shared_queue.front();
			shared_queue.pop();
		}
		pthread_mutex_unlock(&poly_lock);
		// queue is empty but not stopped
		if(!elem){
			usleep(20);
			continue;
		}
		for(MyPolygon *poly:elem->polys){
			ctx->target = poly;
			Pixel *px = poly->getMBB();
			tree.Search(px->low, px->high, MySearchCallback, (void *)ctx);
			if(++local_count==1000){
				pthread_mutex_lock(&report_lock);
				query_count += local_count;
				if(query_count%100000==0){
					log("queried %d polygons",query_count);
				}
				local_count = 0;
				pthread_mutex_unlock(&report_lock);
			}
		}

		if(ctx->gpu){
			elems.push_back(elem);
			if(ctx->candidates.size()>100){
				process_with_gpu(ctx);
				for(queue_element *e:elems){
					delete e;
				}
			}
		}else{
			delete elem;
		}
	}
	if(ctx->gpu){
		process_with_gpu(ctx);
		for(queue_element *e:elems){
			delete e;
		}
	}

	return NULL;
}



int main(int argc, char** argv) {
	string source_path;
	string target_path;
	int num_threads = get_num_threads();
	bool use_partition = false;
	bool use_gpu = false;
	bool in_memory = false;
	int vpr = 10;
	po::options_description desc("query usage");
	desc.add_options()
		("help,h", "produce help message")
		("partition,p", "query with inner partition")
		("gpu,g", "compute with gpu")
		("in_memory,m", "data will be loaded into memory and start")
		("active_partition", "rasterize source polygons before query")
		("source,s", po::value<string>(&source_path), "path to the source")
		("target,t", po::value<string>(&target_path), "path to the target")
		("threads,n", po::value<int>(&num_threads), "number of threads")
		("vpr,v", po::value<int>(&vpr), "number of vertices per raster")
		("element_size,e", po::value<int>(&element_size), "max dimension on vertical")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
		cout << desc << "\n";
		return 0;
	}
	po::notify(vm);
	if(vm.count("partition")){
		use_partition = true;
	}
	if(vm.count("in_memory")){
		in_memory = true;
	}
	if(vm.count("gpu")){
		use_gpu = true;
		gpus = get_gpus();
		for(gpu_info *g:gpus){
			g->init();
		}
		if(gpus.size()==0){
			use_gpu = false;
		}
	}
	if(!vm.count("source")||!vm.count("target")){
		cout << desc << "\n";
		return 0;
	}

	timeval start = get_cur_time();

	vector<MyPolygon *> source = MyPolygon::load_binary_file(source_path.c_str());
	logt("loaded %ld polygons", start, source.size());
	if(use_partition&&vm.count("active_partition")){
		for(MyPolygon *p:source){
			p->partition(vpr);
		}
		logt("partition polygons", start);
	}
	if(use_gpu){
		load_source_togpu(gpus[0], source);
		logt("load to gpu", start);
	}
	for(MyPolygon *p:source){
		tree.Insert(p->getMBB()->low, p->getMBB()->high, p);
	}
	logt("building R-Tree", start);
	pthread_t threads[num_threads];
	query_context ctx[num_threads];
	for(int i=0;i<num_threads;i++){
		ctx[i].thread_id = i;
		ctx[i].vpr = vpr;
		ctx[i].use_partition = use_partition;
		ctx[i].gpu = use_gpu;
	}
	if(!in_memory){
		for(int i=0;i<num_threads;i++){
			pthread_create(&threads[i], NULL, query, (void *)&ctx[i]);
		}
	}

	ifstream is;
	is.open(target_path.c_str(), ios::in | ios::binary);
	int eid = 0;
	queue_element *elem = new queue_element(eid++);
	int loaded = 0;
	size_t pid = 0;
	while(!is.eof()){
		MyPolygon *poly = MyPolygon::read_polygon_binary_file(is);
		if(!poly){
			continue;
		}
		poly->setid(pid++);
		elem->polys.push_back(poly);
		if(elem->polys.size()==element_size){
			while(!in_memory&&shared_queue.size()>10*num_threads){
				usleep(20);
			}
			pthread_mutex_lock(&poly_lock);
			shared_queue.push(elem);
			pthread_mutex_unlock(&poly_lock);
			elem = new queue_element(eid++);
		}
		if(in_memory&&++loaded%100000==0){
			log("loaded %d polygons",loaded);
		}
	}
	if(elem->polys.size()>0){
		pthread_mutex_lock(&poly_lock);
		shared_queue.push(elem);
		pthread_mutex_unlock(&poly_lock);
	}else{
		delete elem;
	}
	stop = true;
	stop = true;

	if(in_memory){
		logt("loaded %d polygons",start,loaded);
		for(int i=0;i<num_threads;i++){
			pthread_create(&threads[i], NULL, query, (void *)&ctx[i]);
		}
	}
	for(int i = 0; i < num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}
	logt("queried %d polygons",start,loaded);

	for(MyPolygon *p:source){
		delete p;
	}
	for(gpu_info *g:gpus){
		delete g;
	}
	return 0;
}



