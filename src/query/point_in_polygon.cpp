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

double *points;
size_t cur_index = 0;
size_t total_points = 0;

size_t inout_checks = 0;
size_t border_checks = 0;
size_t edges_checks = 0;
size_t total_found = 0;

RTree<MyPolygon *, double, 2, double> tree;

bool MySearchCallback(MyPolygon *poly, void* arg){
	query_context *ctx = (query_context *)arg;
	// query with parition
	if(ctx->use_partition){
		if(!poly->ispartitioned()){
			poly->partition(ctx->vpr);
		}
		ctx->found += poly->contain(ctx->target_p,ctx);
	}else if(ctx->use_qtree){
		if(!poly->ispartitioned()){
			 poly->partition_qtree(ctx->vpr);
		}
		QTNode *tnode = poly->get_qtree()->retrieve(ctx->target_p);
		if(tnode->exterior||tnode->interior){
			ctx->found += tnode->interior;
			ctx->inout_check++;
		}else{
			ctx->found += poly->contain(ctx->target_p,ctx);
			ctx->border_check++;
		}
	}else{
		ctx->found += poly->contain(ctx->target_p,ctx);
	}
	return true;
}

int batch_num = 100;
void *query(void *args){
	query_context *ctx = (query_context *)args;
	log("thread %d is started",ctx->thread_id);
	int local_count = 0;
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
			ctx->target_p = Point(points[2*i],points[2*i+1]);
			tree.Search(points+2*i, points+2*i, MySearchCallback, (void *)ctx);
			if(++local_count==1000){
				pthread_mutex_lock(&report_lock);
				query_count += local_count;
				if(query_count%100000==0){
					log("queried %d points",query_count);
				}
				local_count = 0;
				pthread_mutex_unlock(&report_lock);
			}
		}
	}
	pthread_mutex_lock(&report_lock);
	inout_checks += ctx->inout_check;
	border_checks += ctx->border_check;
	edges_checks += ctx->edges_checked;
	total_found += ctx->found;
	pthread_mutex_unlock(&report_lock);

	return NULL;
}



int main(int argc, char** argv) {
	string source_path;
	string target_path;
	int num_threads = get_num_threads();
	int vpr = 10;
	int big_threshold = 500;
	po::options_description desc("query usage");
	desc.add_options()
		("help,h", "produce help message")
		("rasterize,r", "partition with rasterization")
		("qtree,q", "partition with qtree")
		("active_partition,a", "partitioning source polygons before query")
		("source,s", po::value<string>(&source_path), "path to the source")
		("target,t", po::value<string>(&target_path), "path to the target")
		("threads,n", po::value<int>(&num_threads), "number of threads")
		("vpr,v", po::value<int>(&vpr), "number of vertices per raster")
		("big_threshold,b", po::value<int>(&big_threshold), "threshold for complex polygon")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
		cout << desc << "\n";
		return 0;
	}
	po::notify(vm);


	if(!vm.count("source")||!vm.count("target")){
		cout << desc << "\n";
		return 0;
	}

	assert(!(vm.count("rasterize")&&vm.count("qtree")));
	timeval start = get_cur_time();

	vector<MyPolygon *> source = MyPolygon::load_binary_file(source_path.c_str());
	logt("loaded %ld points", start, source.size());

	if(vm.count("rasterize")&&vm.count("active_partition")){
		for(MyPolygon *p:source){
			p->partition(vpr);
		}
		logt("partition polygons", start);
	}

	if(vm.count("qtree")&&vm.count("active_partition")){
		for(MyPolygon *p:source){
			p->partition_qtree(vpr);
		}
		logt("partition polygons", start);
	}

	int treesize = 0;
	for(MyPolygon *p:source){
		if(p->get_num_vertices()>=big_threshold){
			tree.Insert(p->getMBB()->low, p->getMBB()->high, p);
			treesize++;
		}
	}
	logt("building R-Tree with %d nodes", start, treesize);

	// read all the points
	long fsize = file_size(target_path.c_str());
	if(fsize<=0){
		log("%s is empty",target_path.c_str());
		exit(0);
	}else{
		log("size of %s is %ld", target_path.c_str(),fsize);
	}
	total_points = fsize/(2*sizeof(double));

	points = new double[total_points*2];
	ifstream infile(target_path.c_str(), ios::in | ios::binary);
	infile.read((char *)points, fsize);
	logt("loaded %ld points", start,total_points);


	pthread_t threads[num_threads];
	query_context ctx[num_threads];
	for(int i=0;i<num_threads;i++){
		ctx[i].thread_id = i;
		ctx[i].vpr = vpr;
		ctx[i].use_partition = vm.count("rasterize");
		ctx[i].use_qtree = vm.count("qtree");
	}
	for(int i=0;i<num_threads;i++){
		pthread_create(&threads[i], NULL, query, (void *)&ctx[i]);
	}

	for(int i = 0; i < num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}
	logt("queried %d point %ld in/out %ld border %ld edges checked %ld found",start,total_points,inout_checks,border_checks, edges_checks/border_checks, total_found);

	for(MyPolygon *p:source){
		delete p;
	}
	return 0;
}



