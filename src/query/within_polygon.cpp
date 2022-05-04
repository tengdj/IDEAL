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

namespace po = boost::program_options;
using namespace std;


RTree<MyPolygon *, double, 2, double> tree;

bool MySearchCallback(MyPolygon *poly, void* arg){
	query_context *ctx = (query_context *)arg;

	MyPolygon *target= (MyPolygon *)ctx->target;
	// query with rasterization
	if(ctx->use_grid){
		poly->rasterization(ctx->vpr);
		target->rasterization(ctx->vpr);
	}

	if(poly->getid()==target->getid()){
		return true;
	}
	if(poly->getMBB()->distance(*target->getMBB(), ctx->geography)>ctx->distance_buffer_size){
        return true;
	}
	timeval start = get_cur_time();
	ctx->distance = poly->distance(target,ctx);
	ctx->found += ctx->distance<ctx->distance_buffer_size;
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
	pthread_mutex_lock(&gctx->lock);
	log("thread %d is started",ctx->thread_id);
	pthread_mutex_unlock(&gctx->lock);
	ctx->query_count = 0;
	double buffer_low[2];
	double buffer_high[2];

	while(ctx->next_batch(10)){
		for(int i=ctx->index;i<ctx->index_end;i++){
			if(!tryluck(gctx->sample_rate)){
				continue;
			}
			struct timeval query_start = get_cur_time();
			ctx->target = (void *)(gctx->source_polygons[i]);
			box qb = gctx->source_polygons[i]->getMBB()->expand(gctx->distance_buffer_size, true);
			tree.Search(qb.low, qb.high, MySearchCallback, (void *)ctx);
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
    global_ctx.use_grid = true;
    global_ctx.vpr = 20;
    global_ctx.distance_buffer_size = 10;
    global_ctx.geography = true;

	timeval start = get_cur_time();
    global_ctx.source_polygons = MyPolygon::load_binary_file(global_ctx.source_path.c_str(), global_ctx);
	logt("loaded %d objects", start,global_ctx.source_polygons.size());
//
	preprocess(&global_ctx);
	logt("preprocess the source data", start);

	for(MyPolygon *p:global_ctx.source_polygons){
		tree.Insert(p->getMBB()->low, p->getMBB()->high, p);
	}
	logt("building R-Tree with %d nodes", start, global_ctx.source_polygons.size());

	// read all the points
	global_ctx.target_num = 100;
    pthread_t threads[global_ctx.num_threads];
	query_context ctx[global_ctx.num_threads];
	for(int i=0;i<global_ctx.num_threads;i++){
		ctx[i] = global_ctx;//query_context(global_ctx);
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

	return 0;
}



