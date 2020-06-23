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
long query_count = 0;

queue<queue_element *> shared_queue;

RTree<Geometry *, double, 2, double> tree;

int big_threshold = 500;
int small_threshold = 100;
float sample_rate = 1.0;

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
	vector<queue_element *> elems;
	WKTReader *wkt_reader = new WKTReader();
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
			if(sample_rate<1.0 && !tryluck(sample_rate)){
				continue;
			}
			if(poly->get_num_vertices()>small_threshold){
				continue;
			}
			Geometry *geo = wkt_reader->read(poly->to_string());
			ctx->target = (void *)geo;
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
			delete geo;
		}

		delete elem;
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
	bool in_memory = false;
	po::options_description desc("query usage");
	desc.add_options()
		("help,h", "produce help message")
		("in_memory,m", "data will be loaded into memory and start")
		("source,s", po::value<string>(&source_path), "path to the source")
		("target,t", po::value<string>(&target_path), "path to the target")
		("threads,n", po::value<int>(&num_threads), "number of threads")
		("element_size,e", po::value<int>(&element_size), "max dimension on vertical")
		("sample_rate,r", po::value<float>(&sample_rate), "sample rate")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
		cout << desc << "\n";
		return 0;
	}
	po::notify(vm);
	if(vm.count("in_memory")){
		in_memory = true;
	}
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

	vector<MyPolygon *> source = MyPolygon::load_binary_file(source_path.c_str());
	logt("loaded %ld polygons", start, source.size());
	WKTReader *wkt_reader = new WKTReader();
	int treesize = 0;
	for(MyPolygon *p:source){
		if(p->get_num_vertices()>=big_threshold){
			Geometry *geo = wkt_reader->read(p->to_string());
			tree.Insert(p->getMBB()->low, p->getMBB()->high, geo);
			if(++treesize%10000==0){
				log("%d nodes inserted",treesize);
			}
		}
	}
	logt("building R-Tree with %d nodes", start,treesize);
	pthread_t threads[num_threads];
	query_context ctx[num_threads];
	for(int i=0;i<num_threads;i++){
		ctx[i].thread_id = i;
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
	logt("queried %d polygons",start,query_count);

	for(MyPolygon *p:source){
		delete p;
	}
	return 0;
}



