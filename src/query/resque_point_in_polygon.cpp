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
#include "resque_util.h"

namespace po = boost::program_options;
using namespace std;

#define MAX_THREAD_NUM 100

int element_size = 100;

// some shared parameters
pthread_mutex_t poly_lock;
bool stop = false;

pthread_mutex_t report_lock;
long query_count = 0;
long load_count = 0;

double *points;
size_t cur_index = 0;
size_t total_points = 0;
float sample_rate = 1.0;


RTree<Geometry *, double, 2, double> tree;


vector<MyPolygon *> source_polygons;
vector<Geometry *> sources;
vector<Geometry *> targets;

int batch_num = 100;
// multiple thread load
void *load_target(void *args){
	WKTReader *wkt_reader = new WKTReader();
	int local_count = 0;
	char point_buffer[200];
	vector<Geometry *> local_geoms;
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
			sprintf(point_buffer,"POINT(%f %f)",points[2*i],points[2*i+1]);
			Geometry *geo = wkt_reader->read(point_buffer);
			targets[i] = geo;

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
					log("processed %d polygons",load_count);
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

	//exact query
	ctx->found += poly->contains(p);
	ctx->checked++;
	return true;
}

void *query(void *args){
	query_context *ctx = (query_context *)args;
	pthread_mutex_lock(&poly_lock);
	log("thread %d is started",ctx->thread_id);
	pthread_mutex_unlock(&poly_lock);

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
			if(sample_rate<1.0 && !tryluck(sample_rate)){
				continue;
			}
			ctx->target = (void *)(targets[i]);
			ctx->checked = 0;
			ctx->found = 0;
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
	query_count += local_count;
	pthread_mutex_unlock(&report_lock);
	return NULL;
}



int main(int argc, char** argv) {
	string source_path;
	string target_path;
	int num_threads = get_num_threads();
	int big_threshold = 500;
	po::options_description desc("query usage");
	desc.add_options()
		("help,h", "produce help message")
		("source,s", po::value<string>(&source_path)->required(), "path to the source")
		("target,t", po::value<string>(&target_path)->required(), "path to the target")
		("threads,n", po::value<int>(&num_threads), "number of threads")
		("big_threshold,b", po::value<int>(&big_threshold), "threshold for complex polygon")
		("sample_rate,r", po::value<float>(&sample_rate), "sample rate")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
		cout << desc << "\n";
		return 0;
	}
	po::notify(vm);

	timeval start = get_cur_time();
	/////////////////////////////////////////////////////////////////////////////
	// load the source into polygon
	source_polygons = MyPolygon::load_binary_file(source_path.c_str(),big_threshold);
	logt("loaded %ld polygons", start, source_polygons.size());

	/////////////////////////////////////////////////////////////////////////////
	//loading sources as geometry
	sources.resize(source_polygons.size());
	pthread_t threads[num_threads];
	load_count = 0;
	cur_index = 0;
	for(int i=0;i<num_threads;i++){
		pthread_create(&threads[i], NULL, load_source, NULL);
	}
	for(int i = 0; i < num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}
	logt("loading %d polygons into geometry",start,load_count);

	for(int i=0;i<sources.size();i++){
		tree.Insert(source_polygons[i]->getMBB()->low, source_polygons[i]->getMBB()->high, sources[i]);
	}
	logt("building R-Tree with %d nodes", start,sources.size());

	/////////////////////////////////////////////////////////////////////////////////////
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

	////////////////////////////////////////////////////////////////////////////////////
	//loading the points into geometry
	targets.resize(total_points);
	load_count = 0;
	cur_index = 0;
	for(int i=0;i<num_threads;i++){
		pthread_create(&threads[i], NULL, load_target, NULL);
	}

	for(int i = 0; i < num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}
	logt("loading points into geometry",start);

	/////////////////////////////////////////////////////////////////////////////////////
	// querying
	query_context ctx[num_threads];
	load_count = 0;
	cur_index = 0;
	for(int i=0;i<num_threads;i++){
		ctx[i].thread_id = i;
	}
	for(int i=0;i<num_threads;i++){
		pthread_create(&threads[i], NULL, query, (void *)&ctx[i]);
	}

	for(int i = 0; i < num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}
	logt("queried %d points",start,query_count);

	for(MyPolygon *p:source_polygons){
		delete p;
	}
	for(Geometry *g:sources){
		delete g;
	}
	for(Geometry *g:targets){
		delete g;
	}
	delete []points;
	return 0;
}



