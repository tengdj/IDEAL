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

RTree<MyPolygon *, double, 2, double> tree;

bool MySearchCallback(MyPolygon *poly, void* arg){
	query_context *ctx = (query_context *)arg;
	MyPolygon *target = (MyPolygon *)ctx->target;
	// query with parition
	if(ctx->use_grid){
		if(!poly->is_grid_partitioned()){
			poly->partition(ctx->vpr);
		}
		timeval start = get_cur_time();
		bool isfound = poly->contain_try_partition(target->getMBB(), ctx);

		if(ctx->raster_checked_only){
			ctx->pixel_checked++;
			ctx->found += isfound;
		}else{
			ctx->border_checked++;
			ctx->found += poly->contain(target, ctx);
			int nv = poly->get_num_vertices();
			if(nv<5000){
				nv = 100*(nv/100);
				ctx->report_latency(nv, get_time_elapsed(start));
			}
		}
	}else if(ctx->use_qtree){
		if(!poly->is_grid_partitioned()){
			 poly->partition_qtree(ctx->vpr);
		}
		if(poly->get_qtree()->determine_contain(*(target->getMBB()))){
			ctx->pixel_checked++;
		}else{
			ctx->found += poly->contain(target,ctx);
			ctx->border_checked++;
		}
	}else{
		ctx->found += poly->contain(target,ctx);
	}
	// keep going until all hit objects are found
	return true;
}

void *query(void *args){
	query_context *ctx = (query_context *)args;
	query_context *gctx = ctx->global_ctx;
	//log("thread %d is started",ctx->thread_id);
	while(ctx->next_batch(100)){
		for(int i=ctx->index;i<ctx->index_end;i++){
			if(!tryluck(ctx->sample_rate)){
				continue;
			}
			MyPolygon *poly = gctx->target_polygons[i];
			if(poly->get_num_vertices()>100){
				continue;
			}
			ctx->target = (void *)poly;
			Pixel *px = poly->getMBB();
			tree.Search(px->low, px->high, MySearchCallback, (void *)ctx);
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

	global_ctx.source_polygons = MyPolygon::load_binary_file(global_ctx.source_path.c_str(),global_ctx);
	logt("loaded %ld polygons", start, global_ctx.source_polygons.size());

	start = get_cur_time();
	for(MyPolygon *p:global_ctx.source_polygons){
		tree.Insert(p->getMBB()->low, p->getMBB()->high, p);
	}
	logt("building R-Tree with %d nodes", start,global_ctx.source_polygons.size());

	query_context lc(global_ctx);
	lc.big_threshold = 100;
	lc.small_threshold = 0;
	global_ctx.target_polygons = MyPolygon::load_binary_file(global_ctx.target_path.c_str(),lc,true);
	global_ctx.target_num = global_ctx.target_polygons.size();
	logt("loaded %d polygons",start,global_ctx.target_polygons.size());

	int vpr = global_ctx.vpr;
	for(;vpr<=global_ctx.vpr_end;vpr+=10){
		global_ctx.reset_stats();
		global_ctx.vpr = vpr;

		if(global_ctx.use_grid||global_ctx.use_qtree){
			process_partition(&global_ctx);
		}

		pthread_t threads[global_ctx.num_threads];
		query_context ctx[global_ctx.num_threads];
		for(int i=0;i<global_ctx.num_threads;i++){
			ctx[i] = query_context(global_ctx);
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
//		logt("vpr %d: queried %d polygons %ld rastor %ld vector %ld found",start,vpr,global_ctx.query_count,global_ctx.raster_checked,global_ctx.vector_checked
//				,global_ctx.found);
		global_ctx.print_stats();
		logt("query time",start);
		for(MyPolygon *p:global_ctx.source_polygons){
			p->reset_grid_partition();
		}

	}
//	for(MyPolygon *p:source){
//		delete p;
//	}
	return 0;
}



