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
#include "../geos/GEOSTool.h"

namespace po = boost::program_options;
using namespace std;

RTree<Geometry *, double, 2, double> tree;

vector<unique_ptr<Geometry>> sources;
vector<unique_ptr<Geometry>> targets;

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
	query_context *gctx = ctx->global_ctx;
	log("thread %d is started",ctx->thread_id);

	while(ctx->next_batch(100)){
		for(int i=ctx->index;i<ctx->index_end;i++){
			if(!tryluck(ctx->sample_rate)){
				continue;
			}
			ctx->target = (void *)(targets[i].get());
			tree.Search(gctx->target_polygons[i]->getMBB()->low, gctx->target_polygons[i]->getMBB()->high, MySearchCallback, (void *)ctx);
			ctx->report_progress();
		}
	}
	ctx->merge_global();
	return NULL;
}

int main(int argc, char** argv) {
	query_context global_ctx;
	global_ctx = get_parameters(argc, argv);
	timeval start = get_cur_time();
	/////////////////////////////////////////////////////////////////////////////
	// load the source into polygon
	global_ctx.source_polygons = MyPolygon::load_binary_file(global_ctx.source_path.c_str(),global_ctx);
	logt("loaded %ld polygons", start, global_ctx.source_polygons.size());

	/////////////////////////////////////////////////////////////////////////////
	//loading sources as geometry
	process_geometries(&global_ctx,sources);

	for(int i=0;i<sources.size();i++){
		tree.Insert(global_ctx.source_polygons[i]->getMBB()->low, global_ctx.source_polygons[i]->getMBB()->high, sources[i].get());
	}
	logt("building R-Tree with %d nodes", start,sources.size());

	/////////////////////////////////////////////////////////////////////////////////////
	// read all the polygons
	global_ctx.target_polygons = MyPolygon::load_binary_file(global_ctx.target_path.c_str(),global_ctx,true);
	global_ctx.target_num = global_ctx.target_polygons.size();
	////////////////////////////////////////////////////////////////////////////////////
	//loading the polygons into geometry

	process_geometries(&global_ctx,targets,true);

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
	logt("queried %d polygons",start,global_ctx.query_count);

	for(Geometry *g:sources){
		delete g;
	}
	for(Geometry *g:targets){
		delete g;
	}

	return 0;
}
