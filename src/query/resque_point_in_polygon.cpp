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

double *points;
size_t cur_index = 0;
size_t total_points = 0;


RTree<Geometry *, double, 2, double> tree;


bool MySearchCallback(Geometry *poly, void* arg){
	query_context *ctx = (query_context *)arg;
	geos::geom::Geometry *p= (geos::geom::Geometry *)ctx->target;

	//exact query
	ctx->found += poly->covers(p);
	return true;
}

int batch_num = 100;
void *query(void *args){
	query_context *ctx = (query_context *)args;
	pthread_mutex_lock(&poly_lock);
	log("thread %d is started",ctx->thread_id);
	pthread_mutex_unlock(&poly_lock);
	WKTReader *wkt_reader = new WKTReader();
	int local_count = 0;
	char point_buffer[200];
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
			ctx->target = (void *)(geo);
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
			delete geo;
		}
	}
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
		("source,s", po::value<string>(&source_path), "path to the source")
		("target,t", po::value<string>(&target_path), "path to the target")
		("threads,n", po::value<int>(&num_threads), "number of threads")
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

	timeval start = get_cur_time();
//
//	PrecisionModel *pm = new PrecisionModel();
//	GeometryFactory *gf = GeometryFactory::create().get();
	WKTReader *wkt_reader = new WKTReader();

	vector<MyPolygon *> source = MyPolygon::load_binary_file(source_path.c_str());
	logt("loaded %ld points", start, source.size());

	int treesize = 0;
	for(MyPolygon *p:source){
		if(p->get_num_vertices()>=big_threshold){
			Geometry *poly = wkt_reader->read(p->to_string());
			tree.Insert(p->getMBB()->low, p->getMBB()->high, poly);
			treesize++;
		}
	}
	logt("building R-Tree with %d nodes", start,treesize);

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
	}
	for(int i=0;i<num_threads;i++){
		pthread_create(&threads[i], NULL, query, (void *)&ctx[i]);
	}

	for(int i = 0; i < num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}
	logt("queried %d polygons",start,total_points);

	for(MyPolygon *p:source){
		delete p;
	}
	return 0;
}



