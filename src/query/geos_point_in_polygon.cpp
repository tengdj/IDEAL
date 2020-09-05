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

vector<Geometry *> sources;
vector<Geometry *> targets;

bool MySearchCallback(Geometry *poly, void* arg){
	query_context *ctx = (query_context *)arg;
	geos::geom::Geometry *p= (geos::geom::Geometry *)ctx->target;

	assert(poly);
	assert(p);
	//exact query
	ctx->found += poly->contains(p);
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
			ctx->target = (void *)(targets[i]);
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
	timeval start = get_cur_time();
	/////////////////////////////////////////////////////////////////////////////
	// load the source into polygon
	global_ctx.source_polygons = MyPolygon::load_binary_file(global_ctx.source_path.c_str(),global_ctx);

	/////////////////////////////////////////////////////////////////////////////
	//loading sources as geometry
	sources = process_geometries(global_ctx.source_polygons);

	for(int i=0;i<sources.size();i++){
		tree.Insert(global_ctx.source_polygons[i]->getMBB()->low, global_ctx.source_polygons[i]->getMBB()->high, sources[i]);
	}
	logt("building R-Tree with %d nodes", start,sources.size());

	/////////////////////////////////////////////////////////////////////////////////////
	// read all the points
	global_ctx.load_points();

	////////////////////////////////////////////////////////////////////////////////////
	//loading the points into geometry

	targets = process_points(global_ctx.points, global_ctx.target_num);

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
	logt("queried %d points",start,global_ctx.query_count);

	for(Geometry *g:sources){
		delete g;
	}
	for(Geometry *g:targets){
		delete g;
	}

	return 0;
}



