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

double *points;
size_t cur_index = 0;
size_t total_points = 0;
int buffer_size = 10;
float sample_rate = 1.0;

query_context global_ctx;

RTree<MyPolygon *, double, 2, double> tree;
double degree_per_kilometer_latitude = 360.0/40076.0;

double degree_per_kilometer_longitude(double latitude){
	double absla = abs(latitude);
	if(absla==90){
		absla = 89;
	}
	assert(absla<90);
	return 360.0/(sin((90-absla)/90)*40076);
}

bool MySearchCallback(MyPolygon *poly, void* arg){
	query_context *ctx = (query_context *)arg;
	// query with parition
	if(ctx->use_grid){
		if(!poly->is_grid_partitioned()){
			poly->partition(ctx->vpr);
		}
	}
	ctx->distance = poly->distance(ctx->target_p,ctx);
	return true;
}

int batch_num = 100;
void *query(void *args){
	query_context *ctx = (query_context *)args;
	pthread_mutex_lock(&poly_lock);
	log("thread %d is started",ctx->thread_id);
	pthread_mutex_unlock(&poly_lock);
	ctx->query_count = 0;
	double buffer_low[2];
	double buffer_high[2];
	while(cur_index!=total_points){
		int local_cur = 0;
		int local_end = 0;
		pthread_mutex_lock(&poly_lock);
		if(cur_index==total_points){
			pthread_mutex_unlock(&poly_lock);
			break;
		}
		local_cur = cur_index;
		if(local_cur+batch_num>total_points){
			local_end = total_points;
		}else {
			local_end = local_cur+batch_num;
		}
		cur_index = local_end;
		pthread_mutex_unlock(&poly_lock);

		for(int i=local_cur;i<local_end;i++){
			if(sample_rate<1.0&&!tryluck(sample_rate)){
				continue;
			}
			ctx->target_p = Point(points[2*i],points[2*i+1]);
			double shiftx = degree_per_kilometer_longitude(points[2*i+1])*buffer_size;
			double shifty = degree_per_kilometer_latitude*buffer_size;
			buffer_low[0] = points[2*i]-shiftx;
			buffer_low[1] = points[2*i+1]-shifty;
			buffer_high[0] = points[2*i]+shiftx;
			buffer_high[1] = points[2*i+1]+shifty;
			tree.Search(buffer_low, buffer_high, MySearchCallback, (void *)ctx);
			if(++ctx->query_count==1000){
				pthread_mutex_lock(&report_lock);
				global_ctx.query_count += ctx->query_count;
				if(global_ctx.query_count%100000==0){
					log("queried %d points",global_ctx.query_count);
				}
				ctx->query_count = 0;
				pthread_mutex_unlock(&report_lock);
			}
		}
	}
	pthread_mutex_lock(&report_lock);
	global_ctx = global_ctx + *ctx;
	pthread_mutex_unlock(&report_lock);
	return NULL;
}



int main(int argc, char** argv) {
	string source_path;
	string target_path;
	float sample_rate = 1;
	int num_threads = get_num_threads();
	po::options_description desc("query usage");
	desc.add_options()
		("help,h", "produce help message")
		("partition,p", "query with inner partition")
		("source,s", po::value<string>(&source_path), "path to the source")
		("target,t", po::value<string>(&target_path), "path to the target")
		("threads,n", po::value<int>(&num_threads), "number of threads")
		("vpr,v", po::value<int>(&global_ctx.vpr), "number of vertices per raster")
		("buffer_expand,e", po::value<int>(&buffer_size), "buffer in kilometers")
		("big_threshold,b", po::value<int>(&global_ctx.big_threshold), "upper threshold for complex polygon")
		("small_threshold", po::value<int>(&global_ctx.small_threshold), "lower threshold for complex polygon")
		("sample_rate,r", po::value<float>(&sample_rate), "sample rate")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
		cout << desc << "\n";
		return 0;
	}
	po::notify(vm);
	global_ctx.use_grid = vm.count("partition");
	if(!vm.count("source")||!vm.count("target")){
		cout << desc << "\n";
		return 0;
	}

	timeval start = get_cur_time();


	vector<MyPolygon *> source = MyPolygon::load_binary_file(source_path.c_str(), global_ctx);
	logt("loaded %ld points", start, source.size());
	if(global_ctx.use_grid){
		for(MyPolygon *p:source){
			p->partition(global_ctx.vpr);
		}
		logt("partition polygons", start);
	}
	for(MyPolygon *p:source){
		tree.Insert(p->getMBB()->low, p->getMBB()->high, p);
	}
	logt("building R-Tree with %d nodes", start, source.size());

	// read all the points
	long fsize = file_size(target_path.c_str());
	if(fsize<=0){
		log("%s is empty",target_path.c_str());
		exit(0);
	}else{
		log("size of %s is %ld", target_path.c_str(),fsize);
	}
	fsize = fsize*sample_rate;
	fsize = (fsize/(2*sizeof(double)))*2*sizeof(double);
	total_points = fsize/(2*sizeof(double));

	points = new double[total_points*2];
	ifstream infile(target_path.c_str(), ios::in | ios::binary);
	infile.read((char *)points, fsize);
	logt("loaded %ld points", start,total_points);


	pthread_t threads[num_threads];
	query_context ctx[num_threads];
	for(int i=0;i<num_threads;i++){
		ctx[i] = global_ctx;
		ctx[i].thread_id = i;
	}
	for(int i=0;i<num_threads;i++){
		pthread_create(&threads[i], NULL, query, (void *)&ctx[i]);
	}

	for(int i = 0; i < num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}
	logt("queried %d points %d distance calculated %f exterior checked %f boundary checked %d edges checked",start,
			global_ctx.query_count,
			global_ctx.distance_checked,
			1.0*global_ctx.distance_inex_checked/global_ctx.distance_checked,
			1.0*global_ctx.distance_boundary_checked/global_ctx.distance_checked,
			global_ctx.distance_edges_checked/global_ctx.distance_checked);

	for(MyPolygon *p:source){
		delete p;
	}
	return 0;
}



