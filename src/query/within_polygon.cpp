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
#include <boost/program_options.hpp>

namespace po = boost::program_options;
using namespace std;

RTree<Ideal *, double, 2, double> tree;

bool MySearchCallback(Ideal *ideal, void* arg){
	query_context *ctx = (query_context *)arg;

	Ideal *target= (Ideal *)ctx->target;

	if(ideal == target){
		return true;
	}
	// the minimum possible distance is larger than the threshold
	if(ideal->getMBB()->distance(*target->getMBB(), ctx->geography)>ctx->within_distance){
        return true;
	}
	// the maximum possible distance is smaller than the threshold
	if(ideal->getMBB()->max_distance(*target->getMBB(), ctx->geography)<=ctx->within_distance){
		ctx->found++;
        return true;
	}

	ctx->distance = ideal->distance(target,ctx);
	ctx->found += ctx->distance <= ctx->within_distance;

	return true;
}

void *query(void *args){
	query_context *ctx = (query_context *)args;
	query_context *gctx = ctx->global_ctx;
	log("thread %d is started",ctx->thread_id);
	ctx->query_count = 0;
	double buffer_low[2];
	double buffer_high[2];

	while(ctx->next_batch(1)){
		for(int i=ctx->index;i<ctx->index_end;i++){
			if(!tryluck(gctx->sample_rate)){
				ctx->report_progress();
				continue;
			}
			ctx->target = (void *)(gctx->source_ideals[i]);
			box qb = gctx->source_ideals[i]->getMBB()->expand(gctx->within_distance, ctx->geography);
			tree.Search(qb.low, qb.high, MySearchCallback, (void *)ctx);
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
    global_ctx.geography = true;

	timeval start = get_cur_time();
    global_ctx.source_ideals = load_binary_file(global_ctx.source_path.c_str(), global_ctx);
	logt("loaded %d objects", start,global_ctx.source_ideals.size());

	preprocess(&global_ctx);
	logt("preprocess the source data", start);

	for(Ideal *p : global_ctx.source_ideals){
		tree.Insert(p->getMBB()->low, p->getMBB()->high, p);
	}
	logt("building R-Tree with %d nodes", start, global_ctx.source_ideals.size());

	// the target is also the source
	global_ctx.target_num = global_ctx.source_ideals.size();

    pthread_t threads[global_ctx.num_threads];
	query_context ctx[global_ctx.num_threads];
	for(int i=0;i<global_ctx.num_threads;i++){
		ctx[i] = global_ctx; //query_context(global_ctx);
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



