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

bool MySearchCallback(Geometry *poly, void* arg){
	query_context *ctx = (query_context *)arg;
	geos::geom::Geometry *p= (geos::geom::Geometry *)ctx->target;

	ctx->refine_count++;
	ctx->distance = poly->distance(p);
	return true;
}

void *query(void *args){
	query_context *ctx = (query_context *)args;
	query_context *gctx = ctx->global_ctx;
	log("thread %d is started",ctx->thread_id);
	double buffer_low[2];
	double buffer_high[2];
	while(ctx->next_batch(100)){
		for(int i=ctx->index;i<ctx->index_end;i++){
			if(!tryluck(ctx->sample_rate)){
				continue;
			}
			ctx->target = (void *)targets[i].get();
			double shiftx = degree_per_kilometer_longitude(gctx->points[i].y)*gctx->within_distance;
			double shifty = degree_per_kilometer_latitude*gctx->within_distance;
			buffer_low[0] = gctx->points[i].x-shiftx;
			buffer_low[1] = gctx->points[i].y-shifty;
			buffer_high[0] = gctx->points[i].x+shiftx;
			buffer_high[1] = gctx->points[i].y+shifty;
			tree.Search(buffer_low, buffer_high, MySearchCallback, (void *)ctx);
			ctx->report_progress();
		}
	}
	ctx->merge_global();
	return NULL;
}



int main(int argc, char** argv) {
	query_context global_ctx;
	global_ctx = get_parameters(argc,argv);

	WKTReader *wkt_reader = new WKTReader();

	global_ctx.source_polygons = MyPolygon::load_binary_file(global_ctx.source_path.c_str(), global_ctx);

	global_ctx.target_num = global_ctx.source_polygons.size();
	process_geometries(&global_ctx, sources);

	timeval start = get_cur_time();

	for(int i=0;i<sources.size();i++){
		if(sources[i]){
			tree.Insert(global_ctx.source_polygons[i]->getMBB()->low, global_ctx.source_polygons[i]->getMBB()->high, sources[i].get());
		}
	}
	logt("building R-Tree with %d nodes", start,sources.size());

	// read all the points
	global_ctx.load_points();

	process_points(&global_ctx, targets);

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



