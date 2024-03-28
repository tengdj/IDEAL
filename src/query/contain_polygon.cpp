/*
 * Parser.cpp
 *
 *  Created on: May 9, 2020
 *      Author: teng
 */

#include "../include/Ideal.h"
#include <fstream>
#include "../index/RTree.h"
#include <queue>

RTree<Ideal *, double, 2, double> tree;
int ct = 0;

bool MySearchCallback(Ideal *ideal, void* arg){
	query_context *ctx = (query_context *)arg;
	Ideal *target = (Ideal *)ctx->target;
	ctx->found += ideal->contain(target, ctx);
	return true;
}

void *query(void *args){
	query_context *ctx = (query_context *)args;
	query_context *gctx = ctx->global_ctx;
	log("thread %d is started",ctx->thread_id);
	while(ctx->next_batch(10)){
		for(int i=ctx->index;i<ctx->index_end;i++){
			if(!tryluck(ctx->sample_rate)){
				ctx->report_progress();
				continue;
			}
			Ideal *ideal = gctx->target_ideals[i];

			ctx->target = (void *)ideal;
			box *bx = ideal->getMBB();
			tree.Search(bx->low, bx->high, MySearchCallback, (void *)ctx);
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
	global_ctx.source_ideals = load_binary_file(global_ctx.source_path.c_str(),global_ctx);
	start = get_cur_time();
	for(auto p : global_ctx.source_ideals){
		tree.Insert(p->getMBB()->low, p->getMBB()->high, p);
	}
	logt("building R-Tree with %d nodes", start,global_ctx.source_polygons.size());

	global_ctx.target_ideals = load_binary_file(global_ctx.target_path.c_str(),global_ctx);
	global_ctx.target_num = global_ctx.target_ideals.size();
	start = get_cur_time();

	global_ctx.reset_stats();

	preprocess(&global_ctx);
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
//		logt("vpr %d: queried %d polygons %ld rastor %ld vector %ld found",start,vpr,global_ctx.query_count,global_ctx.raster_checked,global_ctx.vector_checked
//				,global_ctx.found);
	global_ctx.print_stats();
	logt("query",start);

//	for(MyPolygon *p:source){
//		delete p;
//	}
	return 0;
}



