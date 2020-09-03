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

#include <geos/geom/PrecisionModel.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h>
#include <geos/operation/distance/DistanceOp.h>
#include <geos/geom/Point.h>
#include <geos/io/WKTReader.h>
#include <geos/io/WKTWriter.h>
#include <geos/opBuffer.h>

using namespace geos;
using namespace geos::io;
using namespace geos::geom;
using namespace geos::operation::buffer;
using namespace geos::operation::distance;

namespace po = boost::program_options;
using namespace std;

#define MAX_THREAD_NUM 100

int element_size = 100;

// some shared parameters
pthread_mutex_t poly_lock;
bool stop = false;

pthread_mutex_t report_lock;
queue<queue_element *> shared_queue;

RTree<Geometry *, double, 2, double> tree;

vector<MyPolygon *> source_polygons;
vector<MyPolygon *> target_polygons;

vector<Geometry *> sources;
vector<Geometry *> targets;

int cur_index = 0;
int load_count = 0;
query_context global_ctx;


int batch_num = 100;
// multiple thread load
void *load_target(void *args){
	WKTReader *wkt_reader = new WKTReader();
	int local_count = 0;
	vector<Geometry *> local_geoms;
	while(cur_index!=target_polygons.size()){
		int local_cur = 0;
		int local_end = 0;
		pthread_mutex_lock(&poly_lock);
		if(cur_index==target_polygons.size()){
			pthread_mutex_unlock(&poly_lock);
			break;
		}
		local_cur = cur_index;
		if(local_cur+batch_num>target_polygons.size()){
			local_end = target_polygons.size();
		}else {
			local_end = local_cur+batch_num;
		}
		cur_index = local_end;
		pthread_mutex_unlock(&poly_lock);

		for(int i=local_cur;i<local_end;i++){
			targets[i] = wkt_reader->read(target_polygons[i]->to_string());
			if(++local_count==1000){
				pthread_mutex_lock(&report_lock);
				load_count += local_count;
				if(load_count%1000000==0){
					log("loaded %d points",load_count);
				}
				local_count = 0;
				pthread_mutex_unlock(&report_lock);
			}
		}
	}
	pthread_mutex_lock(&report_lock);
	load_count += local_count;
	pthread_mutex_unlock(&report_lock);
	delete wkt_reader;
	return NULL;
}

// multiple thread load
void *load_source(void *args){
	WKTReader *wkt_reader = new WKTReader();
	int local_count = 0;
	while(cur_index!=source_polygons.size()){
		int local_cur = 0;
		int local_end = 0;
		pthread_mutex_lock(&poly_lock);
		if(cur_index==source_polygons.size()){
			pthread_mutex_unlock(&poly_lock);
			break;
		}
		local_cur = cur_index;
		if(local_cur+batch_num>source_polygons.size()){
			local_end = source_polygons.size();
		}else {
			local_end = local_cur+batch_num;
		}
		cur_index = local_end;
		pthread_mutex_unlock(&poly_lock);

		for(int i=local_cur;i<local_end;i++){
			Geometry *geo = wkt_reader->read(source_polygons[i]->to_string());
			sources[i] = geo;
			if(++local_count==1000){
				pthread_mutex_lock(&report_lock);
				load_count += local_count;
				if(load_count%100000==0){
					log("load %d polygons",load_count);
				}
				local_count = 0;
				pthread_mutex_unlock(&report_lock);
			}
		}
	}
	pthread_mutex_lock(&report_lock);
	load_count += local_count;
	pthread_mutex_unlock(&report_lock);
	delete wkt_reader;
	return NULL;
}

bool MySearchCallback(Geometry *poly, void* arg){
	query_context *ctx = (query_context *)arg;
	geos::geom::Geometry *p= (geos::geom::Geometry *)ctx->target;
	try{
		ctx->found += poly->covers(p);
	}catch(...){

	}
	// keep going until all hit objects are found
	return true;
}

void *query(void *args){
	query_context *ctx = (query_context *)args;
	pthread_mutex_lock(&poly_lock);
	log("thread %d is started",ctx->thread_id);
	pthread_mutex_unlock(&poly_lock);

	int local_count = 0;
	while(cur_index!=target_polygons.size()){
		int local_cur = 0;
		int local_end = 0;
		pthread_mutex_lock(&poly_lock);
		if(cur_index==target_polygons.size()){
			pthread_mutex_unlock(&poly_lock);
			break;
		}
		local_cur = cur_index;
		if(local_cur+batch_num>target_polygons.size()){
			local_end = target_polygons.size();
		}else {
			local_end = local_cur+batch_num;
		}
		cur_index = local_end;
		pthread_mutex_unlock(&poly_lock);

		for(int i=local_cur;i<local_end;i++){
			ctx->target = (void *)(targets[i]);
			tree.Search(target_polygons[i]->getMBB()->low, target_polygons[i]->getMBB()->high, MySearchCallback, (void *)ctx);
			if(++local_count==1000){
				pthread_mutex_lock(&report_lock);
				global_ctx.query_count += local_count;
				if(global_ctx.query_count%100000==0){
					log("queried %d points",global_ctx.query_count);
				}
				local_count = 0;
				pthread_mutex_unlock(&report_lock);
			}
		}
	}
	ctx->merge_global();
	return NULL;
}



int main(int argc, char** argv) {
	string source_path;
	string target_path;
	int num_threads = get_num_threads();
	po::options_description desc("query usage");
	desc.add_options()
		("help,h", "produce help message")
		("source,s", po::value<string>(&source_path), "path to the source")
		("target,t", po::value<string>(&target_path), "path to the target")
		("threads,n", po::value<int>(&num_threads), "number of threads")
		("big_threshold,b", po::value<int>(&global_ctx.big_threshold), "up threshold for complex polygon")
		("small_threshold", po::value<int>(&global_ctx.small_threshold), "low threshold for complex polygon")
		("sample_rate,r", po::value<float>(&global_ctx.sample_rate), "sample rate")
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
	if(!file_exist(source_path.c_str())){
		log("%s does not exist",source_path.c_str());
		return 0;
	}
	if(!file_exist(target_path.c_str())){
		log("%s does not exist",target_path.c_str());
		return 0;
	}
	timeval start = get_cur_time();
	/*---------------------------------------------------------------------*/
	global_ctx.sort_polygons = true;
	source_polygons = MyPolygon::load_binary_file(source_path.c_str(),global_ctx);
	logt("loaded %ld polygons", start, source_polygons.size());
	/*---------------------------------------------------------------------*/
	sources.resize(source_polygons.size());
	pthread_t threads[num_threads];
	query_context ctx[num_threads];
	for(int i=0;i<num_threads;i++){
		ctx[i] = global_ctx;
		ctx[i].thread_id = i;
	}

	load_count = 0;
	cur_index = 0;
	for(int i=0;i<num_threads;i++){
		pthread_create(&threads[i], NULL, load_source, NULL);
	}
	for(int i = 0; i < num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}
	logt("parsed %d source polygons into geometry",start,load_count);

	for(int i=0;i<sources.size();i++){
		tree.Insert(source_polygons[i]->getMBB()->low, source_polygons[i]->getMBB()->high, sources[i]);
	}
	logt("built R-Tree with %d nodes", start,sources.size());
	/*---------------------------------------------------------------------*/

	query_context lc = global_ctx;
	lc.big_threshold = 100;
	lc.small_threshold = 0;
	target_polygons = MyPolygon::load_binary_file(target_path.c_str(),lc);
	logt("loaded %d polygons",start,target_polygons.size());
	load_count = 0;
	cur_index = 0;
	targets.resize(target_polygons.size());
	for(int i=0;i<num_threads;i++){
		pthread_create(&threads[i], NULL, load_target, NULL);
	}
	for(int i = 0; i < num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}
	logt("parsed %d target polygons into geometry",start,load_count);
	/*---------------------------------------------------------------------*/
	cur_index = 0;
	for(int i=0;i<num_threads;i++){
		pthread_create(&threads[i], NULL, query, (void *)&ctx[i]);
	}

	for(int i = 0; i < num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}
	logt("queried %d polygons",start,global_ctx.found);

	for(MyPolygon *p:source_polygons){
		delete p;
	}
	for(MyPolygon *p:target_polygons){
		delete p;
	}
	return 0;
}



