/*
 * Parser.cpp
 *
 *  Created on: May 9, 2020
 *      Author: teng
 */



#include "../index/RTree.h"
#include <queue>
#include <fstream>
#include "../geometry/MyPolygon.h"



// some shared parameters

RTree<MyPolygon *, double, 2, double> tree;

bool MySearchCallback(MyPolygon *poly, void* arg){
	query_context *ctx = (query_context *)arg;

	Point p = *(Point *)ctx->target;
	// query with parition
	if(ctx->use_grid){
		assert(poly->is_grid_partitioned());
		struct timeval start = get_cur_time();
		ctx->found += poly->contain(p,ctx);
		if(!ctx->raster_checked_only){
			int nv = poly->get_num_vertices();
			if(nv<5000){
				nv = 100*(nv/100);
				ctx->report_latency(nv, get_time_elapsed(start));
			}
		}
	}else if(ctx->use_qtree){
		assert(poly->is_qtree_partitioned());
		QTNode *tnode = poly->get_qtree()->retrieve(p);
		if(tnode->exterior||tnode->interior){
			ctx->found += tnode->interior;
			ctx->pixel_checked++;
		}else if(ctx->query_vector){
			ctx->found += poly->contain(p,ctx);
			ctx->border_checked++;
		}
	}else{
		ctx->found += poly->contain(p,ctx);
		ctx->border_checked++;
	}
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
			Point p(gctx->points[2*i],gctx->points[2*i+1]);
			ctx->target = (void *)&p;
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
	global_ctx.query_type = QueryType::contain;

	global_ctx.sort_polygons = true;
	global_ctx.source_polygons = MyPolygon::load_binary_file(global_ctx.source_path.c_str(),global_ctx);

	timeval start = get_cur_time();
	if(global_ctx.use_grid||global_ctx.use_qtree){
		process_partition(&global_ctx);
		start = get_cur_time();
		if(!global_ctx.query_vector){
			//collect partitioning status
			size_t total_data_size = 0;
			size_t total_partition_size = 0;
			for(MyPolygon *poly:global_ctx.source_polygons){
				total_data_size += poly->get_data_size();
				if(global_ctx.use_qtree){
					total_partition_size += poly->get_qtree()->size();
				}
			}
			log("%ld %ld",total_data_size,total_partition_size);
			return 0;
		}
	}

	start = get_cur_time();
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



