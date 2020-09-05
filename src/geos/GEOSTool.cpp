/*
 * GEOSTool.cpp
 *
 *  Created on: Sep 4, 2020
 *      Author: teng
 */


#include "GEOSTool.h"



void *load_source(void *args){
	query_context *ctx = (query_context *)args;
	query_context *gctx = ctx->global_ctx;
	log("thread %d started",ctx->thread_id);
	vector<Geometry *> *dest = (vector<Geometry *> *)gctx->target;
	WKTReader *wkt_reader = new WKTReader();
	while(ctx->next_batch(100)){
		for(int i=ctx->index;i<ctx->index_end;i++){
			(*dest)[i] = wkt_reader->read(gctx->source_polygons[i]->to_string());
		}
	}
	delete wkt_reader;
	return NULL;
}


vector<Geometry *> process_geometries(vector<MyPolygon *> &polys){
	struct timeval start = get_cur_time();
	vector<Geometry *> dest;
	dest.resize(polys.size());
	query_context global_ctx;
	global_ctx.target = &dest;
	global_ctx.source_polygons = polys;
	global_ctx.target_num = polys.size();

	pthread_t threads[global_ctx.num_threads];
	query_context ctx[global_ctx.num_threads];
	for(int i=0;i<global_ctx.num_threads;i++){
		ctx[i] = global_ctx;
		ctx[i].thread_id = i;
		ctx[i].global_ctx = &global_ctx;
	}

	for(int i=0;i<global_ctx.num_threads;i++){
		pthread_create(&threads[i], NULL, load_source, (void *)&ctx[i]);
	}

	for(int i = 0; i < global_ctx.num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}

	logt("%d polygons are processed",start,dest.size());
	return dest;
}


void *load_points(void *args){
	query_context *ctx = (query_context *)args;
	query_context *gctx = ctx->global_ctx;

	vector<Geometry *> *dest = (vector<Geometry *> *)gctx->target;
	WKTReader *wkt_reader = new WKTReader();
	char point_buffer[200];
	while(ctx->next_batch(100)){
		for(int i=ctx->index;i<ctx->index_end;i++){
			sprintf(point_buffer,"POINT(%f %f)",gctx->points[2*i],gctx->points[2*i+1]);
			(*dest)[i] = wkt_reader->read(point_buffer);
		}
	}
	delete wkt_reader;
	return NULL;
}


vector<Geometry *> process_points(double *points, size_t point_num){

	struct timeval start = get_cur_time();
	vector<Geometry *> ret;
	ret.resize(point_num);
	query_context global_ctx;
	global_ctx.target = &ret;
	global_ctx.points = points;
	global_ctx.target_num = point_num;

	pthread_t threads[global_ctx.num_threads];
	query_context ctx[global_ctx.num_threads];
	for(int i=0;i<global_ctx.num_threads;i++){
		ctx[i] = global_ctx;
		ctx[i].thread_id = i;
		ctx[i].global_ctx = &global_ctx;
	}

	for(int i=0;i<global_ctx.num_threads;i++){
		pthread_create(&threads[i], NULL, load_points, (void *)&ctx[i]);
	}

	for(int i = 0; i < global_ctx.num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}
	logt("%d points are processed",start,ret.size());
	return ret;
}


