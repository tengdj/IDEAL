/*
 * Parser.cpp
 *
 *  Created on: May 9, 2020
 *      Author: teng
 */



#include "../index/RTree.h"
#include <queue>
#include <fstream>
#include "../include/Ideal.h"



// some shared parameters

RTree<Ideal *, double, 2, double> ideal_rtree;
RTree<MyPolygon *, double, 2, double> poly_rtree;

bool IdealSearchCallback(Ideal *ideal, void* arg){
	query_context *ctx = (query_context *)arg;
	ctx->found += ideal->contain(*(Point *)ctx->target, ctx);
	return true;
}

bool PolygonSearchCallback(MyPolygon *poly, void* arg){
	query_context *ctx = (query_context *)arg;
	ctx->found += poly->contain(*(Point *)ctx->target, ctx);
	return true;
}

void *query(void *args){
	query_context *ctx = (query_context *)args;
	query_context *gctx = ctx->global_ctx;
	log("thread %d is started",ctx->thread_id);
	char point_buffer[200];

	while(ctx->next_batch(100)){
		for(int i=ctx->index;i<ctx->index_end;i++){
			if(!tryluck(ctx->sample_rate)){
				ctx->report_progress();
				continue;
			}
			ctx->target = (void *)&gctx->points[i];
			if(gctx->use_ideal){
				ideal_rtree.Search((double *)(gctx->points+i), (double *)(gctx->points+i), IdealSearchCallback, (void *)ctx);
			}
			if(gctx->use_vector){
				poly_rtree.Search((double *)(gctx->points+i), (double *)(gctx->points+i), PolygonSearchCallback, (void *)ctx);
			}
			ctx->report_progress();
		}
	}
	ctx->merge_global();

	return NULL;
}



int main(int argc, char** argv) {

	query_context global_ctx;
	global_ctx = get_parameters(argc, argv);
	global_ctx.query_type = QueryType::contain;

	if(global_ctx.use_ideal){
		global_ctx.source_ideals = load_binary_file(global_ctx.source_path.c_str(),global_ctx);
		timeval start = get_cur_time();
		for(auto p : global_ctx.source_ideals){
			ideal_rtree.Insert(p->getMBB()->low, p->getMBB()->high, p);
		}
		logt("building R-Tree with %d nodes", start, global_ctx.source_ideals.size());
	}else{
		global_ctx.source_polygons = load_polygons_from_path(global_ctx.source_path.c_str(),global_ctx);
		timeval start = get_cur_time();
		for(auto p : global_ctx.source_polygons){
			poly_rtree.Insert(p->getMBB()->low, p->getMBB()->high, p);
		}
		logt("building R-Tree with %d nodes", start, global_ctx.source_polygons.size());
	}

	preprocess(&global_ctx);

	// read all the points
	global_ctx.load_points();


	timeval start = get_cur_time();
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


	return 0;
}



