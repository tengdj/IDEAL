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
	// query with parition
	if(ctx->use_grid){
		if(!poly->is_grid_partitioned()){
			poly->partition(ctx->vpr);
		}
	}
	Point *p = (Point *)ctx->target;
	if(poly->getMBB()->distance_geography(*p)>(double)ctx->distance_buffer_size){
        return true;
	}

	timeval start = get_cur_time();
	ctx->distance = poly->distance(*p,ctx);
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

	while(ctx->next_batch(100)){

		for(int i=ctx->index;i<ctx->index_end;i++){
			if(!tryluck(gctx->sample_rate)){
				continue;
			}
			struct timeval query_start = get_cur_time();
			Point p(gctx->points[2*i],gctx->points[2*i+1]);
			ctx->target = (void *)&p;
			double shiftx = degree_per_kilometer_longitude(gctx->points[2*i+1])*gctx->distance_buffer_size;
			double shifty = degree_per_kilometer_latitude*gctx->distance_buffer_size;
			buffer_low[0] = gctx->points[2*i]-shiftx;
			buffer_low[1] = gctx->points[2*i+1]-shifty;
			buffer_high[0] = gctx->points[2*i]+shiftx;
			buffer_high[1] = gctx->points[2*i+1]+shifty;
			tree.Search(buffer_low, buffer_high, MySearchCallback, (void *)ctx);
			ctx->report_progress();
			ctx->check_time += get_time_elapsed(query_start);
		}
	}
	ctx->merge_global();
	return NULL;
}



int main(int argc, char** argv) {
	query_context global_ctx;
	global_ctx = get_parameters(argc, argv);
	global_ctx.query_type = QueryType::within;

	timeval start = get_cur_time();

    global_ctx.source_polygons = MyPolygon::load_binary_file(global_ctx.source_path.c_str(), global_ctx);
	if(global_ctx.use_grid||global_ctx.use_qtree){
		process_partition(&global_ctx);
	}

	if(global_ctx.use_convex_hull){
		process_convex_hull(&global_ctx);
	}

	if(global_ctx.use_mer){
		process_mer(&global_ctx);
	}


	if(global_ctx.use_triangulate){
		if(global_ctx.valid_path.size()>0){
			 ifstream is(global_ctx.valid_path);
			 int num = 0;
			 while(is>>num){
				 assert(num<global_ctx.source_polygons.size());
				 global_ctx.source_polygons[num]->valid_for_triangulate = true;
			 }
			 is.close();
		}
		process_triangulate(&global_ctx);
		process_internal_rtree(&global_ctx);
	}

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



