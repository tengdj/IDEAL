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

namespace po = boost::program_options;
using namespace std;


RTree<MyPolygon *, double, 2, double> tree;

bool MySearchCallback(MyPolygon *poly, void* arg){
	query_context *ctx = (query_context *)arg;
	// query with parition
	if(ctx->use_grid){
		poly->rasterization(ctx->vpr);
	}
	Point *p = (Point *)ctx->target;
	if(poly->getMBB()->distance(*p, ctx->geography)>ctx->within_distance){
        return true;
	}

	timeval start = get_cur_time();
	ctx->distance = poly->distance(*p,ctx);
	ctx->found += ctx->distance <= ctx->within_distance;
	if(ctx->collect_latency){
		int nv = poly->get_num_vertices();
		if(nv<5000){
			nv = 100*(nv/100);
			ctx->report_latency(nv, get_time_elapsed(start));
		}
	}
	return true;
}

void *query(void *args){
	query_context *ctx = (query_context *)args;
	query_context *gctx = ctx->global_ctx;
	gctx->lock();
	log("thread %d is started",ctx->thread_id);
	gctx->unlock();
	ctx->query_count = 0;
	double buffer_low[2];
	double buffer_high[2];

	while(ctx->next_batch(100)){
		for(int i=ctx->index;i<ctx->index_end;i++){
			if(!tryluck(gctx->sample_rate)){
				ctx->report_progress();
				continue;
			}
			struct timeval query_start = get_cur_time();
			ctx->target = (void *)&gctx->points[i];
			double shiftx = degree_per_kilometer_longitude(gctx->points[i].y)*gctx->within_distance;
			double shifty = degree_per_kilometer_latitude*gctx->within_distance;
			buffer_low[0] = gctx->points[i].x-shiftx;
			buffer_low[1] = gctx->points[i].y-shifty;
			buffer_high[0] = gctx->points[i].x+shiftx;
			buffer_high[1] = gctx->points[i].y+shifty;
			tree.Search(buffer_low, buffer_high, MySearchCallback, (void *)ctx);
			ctx->object_checked.execution_time += get_time_elapsed(query_start);
			ctx->report_progress();
		}
	}
	ctx->merge_global();
	return NULL;
}



int main(int argc, char** argv) {
	query_context global_ctx;
	global_ctx = get_parameters(argc, argv);
	global_ctx.query_type = QueryType::within;


    global_ctx.source_polygons = load_binary_file(global_ctx.source_path.c_str(), global_ctx);

	preprocess(&global_ctx);

	timeval start = get_cur_time();
	for(MyPolygon *p:global_ctx.source_polygons){
		tree.Insert(p->getMBB()->low, p->getMBB()->high, p);
	}
	logt("building R-Tree with %d nodes", start, global_ctx.source_polygons.size());

	// read all the points
	global_ctx.load_points();

	start = get_cur_time();

    pthread_t threads[global_ctx.num_threads];
	query_context ctx[global_ctx.num_threads];
	for(int i=0;i<global_ctx.num_threads;i++){
		ctx[i] = query_context(global_ctx);
		ctx[i].thread_id = i;
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

	return 0;
}



