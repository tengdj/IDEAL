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

namespace po = boost::program_options;
using namespace std;

#define MAX_THREAD_NUM 100

int element_size = 100;

// some shared parameters
pthread_mutex_t poly_lock;
bool stop = false;

pthread_mutex_t report_lock;
long query_count = 0;

queue<queue_element *> shared_queue;

RTree<MyPolygon *, double, 2, double> tree;


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
	}
	ctx->found += poly->contain(ctx->target, ctx);
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

		delete elem;
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
	int count100 = 0;
	int count500 = 0;
	while(!is.eof()){
		MyPolygon *poly = MyPolygon::read_polygon_binary_file(is);
		if(!poly){
			continue;
		}
		if(poly->get_num_vertices()>100){
			count100++;
		}
		if(poly->get_num_vertices()>500){
			count500++;
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

	cout<<count100<<" "<<count500<<endl;
	for(MyPolygon *p:source){
		delete p;
	}
	return 0;
}



