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
	Point p = *(Point *)ctx->target;
	if(poly->getMBB()->distance_geography(p)>ctx->distance_buffer_size){
        //printf("%f\n",poly->getMBB()->distance_geography(p));
        return true;
	}

	timeval start = get_cur_time();
	ctx->distance = poly->distance(p,ctx);
	int nv = poly->get_num_vertices();
	if(nv<5000){
		nv = 100*(nv/100);
		ctx->report_latency(nv, get_time_elapsed(start));
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
		}
	}
	ctx->merge_global();
	return NULL;
}



int main(int argc, char** argv) {
	query_context global_ctx;
	global_ctx = get_parameters(argc, argv);
	global_ctx.query_type = QueryType::within;
	bool use_qtree = global_ctx.use_qtree;

	timeval start = get_cur_time();

    global_ctx.source_polygons = MyPolygon::load_binary_file(global_ctx.source_path.c_str(), global_ctx);
	logt("loaded %ld polygons", start, global_ctx.source_polygons.size());
	if(use_qtree){
        global_ctx.use_grid = true;
        global_ctx.use_qtree = false;
	}
	if(global_ctx.use_grid){
		process_partition(&global_ctx);
	}
    global_ctx.use_qtree = use_qtree;

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
	logt("queried %ld points %ld distance calculated %f exterior checked %f boundary checked %ld edges checked %ld edges per boundary pixel",start,
			global_ctx.query_count,
			global_ctx.checked_count,
			1.0*global_ctx.raster_checked/global_ctx.checked_count,
			1.0*global_ctx.vector_checked/global_ctx.checked_count,
			global_ctx.edges_checked/global_ctx.checked_count,
			global_ctx.edges_checked/global_ctx.vector_checked);
	for(auto it:global_ctx.vertex_number){
		cout<<it.first<<"\t"<<global_ctx.latency[it.first]/it.second<<endl;
	}
	return 0;
}



