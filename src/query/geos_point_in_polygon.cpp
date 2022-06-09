/*
 * Parser.cpp
 *
 *  Created on: May 9, 2020
 *      Author: teng
 */

#include "../include/MyPolygon.h"
#include <fstream>
#include "../index/RTree.h"
#include <queue>
#include <boost/program_options.hpp>
#include "../include/GEOSTool.h"

namespace po = boost::program_options;
using namespace std;

RTree<Geometry *, double, 2, double> tree;

vector<unique_ptr<Geometry>> sources;
vector<unique_ptr<Geometry>> targets;

int idx = 0;
bool MySearchCallback(Geometry *poly, void* arg){
	query_context *ctx = (query_context *)arg;
	geos::geom::Geometry *p= (geos::geom::Geometry *)ctx->target;

	assert(poly);
	assert(p);
	//exact query
	ctx->refine_count++;
	//bool found = (geos::geom::Polygon *)poly->contains((geos::geom::Point *)p)==0;
	bool found = poly->distance(p)==0;

	ctx->found += found;
//	double dist = poly->distance(p);
//	if(found){
//		assert(dist==0);
//	}else{
//		if(dist>0){
//			cout<<poly->toString()<<endl;
//			cout<<p->toString()<<endl;
//			cout<<poly->getArea()<<endl;
//			printf("%f\n",dist);
//			assert(0);
//		}
//	}

	return true;
}

void *query(void *args){
	query_context *ctx = (query_context *)args;
	query_context *gctx = ctx->global_ctx;
	log("thread %d is started",ctx->thread_id);

	while(ctx->next_batch(100)){
		for(int i=ctx->index;i<ctx->index_end;i++){
			if(!tryluck(ctx->sample_rate)){
				continue;
			}
			ctx->target = (void *)(targets[i].get());
			tree.Search(gctx->points+2*i, gctx->points+2*i, MySearchCallback, (void *)ctx);
			ctx->report_progress();
		}
	}
	ctx->merge_global();
	return NULL;
}



int main(int argc, char** argv) {
	query_context global_ctx;
	global_ctx = get_parameters(argc, argv);
	/////////////////////////////////////////////////////////////////////////////
	// load the source into polygon
	global_ctx.source_polygons = MyPolygon::load_binary_file(global_ctx.source_path.c_str(),global_ctx);

	/////////////////////////////////////////////////////////////////////////////
	//loading sources as geometry
	global_ctx.target_num = global_ctx.source_polygons.size();
	process_geometries(&global_ctx, sources);

	timeval start = get_cur_time();
	int idx = 0;
	for(int i=0;i<sources.size();i++){
		if(sources[i]){
			assert(sources[i].get()->isPolygonal());

			box *mbr = global_ctx.source_polygons[i]->getMBB();
			tree.Insert(mbr->low, mbr->high, sources[i].get());
			idx++;
		}
	}
	logt("building R-Tree with %d nodes", start, idx);

	/////////////////////////////////////////////////////////////////////////////////////
	// read all the points
	global_ctx.load_points();

	////////////////////////////////////////////////////////////////////////////////////
	//loading the points into geometry

	process_points(&global_ctx,targets);

	/////////////////////////////////////////////////////////////////////////////////////
	// querying
	start = get_cur_time();
	pthread_t threads[global_ctx.num_threads];
	query_context ctx[global_ctx.num_threads];
	for(int i=0;i<global_ctx.num_threads;i++){
		ctx[i] = global_ctx;
		ctx[i].thread_id = i;
		ctx[i].global_ctx = &global_ctx;
	}
	for(int i=0;i<global_ctx.num_threads;i++){
		pthread_create(&threads[i], NULL, query, (void *)&ctx[i]);
	}

	for(int i = 0; i < global_ctx.num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}
	global_ctx.print_stats();
	logt("total query",start);
	sources.clear();
	targets.clear();
	return 0;
}



