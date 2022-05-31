/*
 * GEOSTool.cpp
 *
 *  Created on: Sep 4, 2020
 *      Author: teng
 */


#include "GEOSTool.h"


void *load_source(void *args){
	query_context *ctx = (query_context *)args;
	log("thread %d started to load GEOS object",ctx->thread_id);
	vector<unique_ptr<Geometry>> *dest = (vector<unique_ptr<Geometry>> *)ctx->target;
	WKTReader *wkt_reader = new WKTReader();
	while(ctx->next_batch(10)){
		for(int i=ctx->index;i<ctx->index_end;i++){
			try{
				//log("processing %d", ctx->source_polygons[i]->getid());
				(*dest)[i] = wkt_reader->read(ctx->source_polygons[i]->to_string(false, true));
				(*dest)[i].get()->normalize();
				//log("%f",(*dest)[i]->getArea());
			}catch(...){
				log("failed to parse polygon %ld",ctx->source_polygons[i]->getid());
				(*dest)[i] = NULL;
			}
			ctx->report_progress();
		}
	}
	delete wkt_reader;
	return NULL;
}


void process_geometries(query_context *global_ctx, vector<unique_ptr<Geometry>> &dest, bool process_target){
	struct timeval start = get_cur_time();

	vector<MyPolygon *> polys = process_target?global_ctx->target_polygons:global_ctx->source_polygons;
	dest.resize(polys.size());

	int num_threads = global_ctx->num_threads;
	global_ctx->index = 0;
	global_ctx->target_num = polys.size();
	pthread_t threads[num_threads];
	query_context ctx[num_threads];
	for(int i=0;i<num_threads;i++){
		ctx[i] = *global_ctx;
		ctx[i].thread_id = i;
		ctx[i].global_ctx = global_ctx;
		ctx[i].target = &dest;
		ctx[i].source_polygons = polys;
		ctx[i].target_num = polys.size();
	}

	for(int i=0;i<num_threads;i++){
		pthread_create(&threads[i], NULL, load_source, (void *)&ctx[i]);
	}

	for(int i = 0; i < num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}
	global_ctx->reset_stats();

	logt("%d polygons are processed",start,dest.size());
}


void *load_points(void *args){
	query_context *ctx = (query_context *)args;

	vector<unique_ptr<Geometry>> *dest = (vector<unique_ptr<Geometry>> *)ctx->target;
	WKTReader *wkt_reader = new WKTReader();
	char point_buffer[200];
	while(ctx->next_batch(100)){
		for(int i=ctx->index;i<ctx->index_end;i++){
			sprintf(point_buffer,"POINT(%f %f)",ctx->points[2*i],ctx->points[2*i+1]);
			(*dest)[i] = wkt_reader->read(point_buffer);
		}
	}
	delete wkt_reader;
	return NULL;
}


void process_points(query_context *global_ctx,vector<unique_ptr<Geometry>> &dest){

	struct timeval start = get_cur_time();
	dest.resize(global_ctx->target_num);
	int num_threads = global_ctx->num_threads;

	pthread_t threads[num_threads];
	query_context ctx[num_threads];
	for(int i=0;i<num_threads;i++){
		ctx[i] = *global_ctx;
		ctx[i].thread_id = i;
		ctx[i].global_ctx = global_ctx;
		ctx[i].target = &dest;
		ctx[i].points = global_ctx->points;
		ctx[i].target_num = global_ctx->target_num;
	}

	for(int i=0;i<num_threads;i++){
		pthread_create(&threads[i], NULL, load_points, (void *)&ctx[i]);
	}

	for(int i = 0; i < num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}
	global_ctx->index = 0;
	logt("%d points are processed",start,dest.size());
}


