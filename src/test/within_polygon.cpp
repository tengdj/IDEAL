/*
 * Parser.cpp
 *
 *  Created on: May 9, 2020
 *      Author: teng
 */

#include "../include/MyPolygon.h"
#include <fstream>
#include "../index/RTree.h"
#include <queue>
#include <boost/program_options.hpp>

namespace po = boost::program_options;
using namespace std;

RTree<MyPolygon *, double, 2, double> tree;

size_t total_checked = 0;

size_t judge_count[5] = {0,0,0,0,0};
double right_total[5] = {0,0,0,0,0};

double wrong_total[5] = {0,0,0,0,0};

double optimal = 0.0;
pthread_mutex_t judge_t;

bool MySearchCallback(MyPolygon *poly, void* arg){
	query_context *ctx = (query_context *)arg;

	MyPolygon *target= (MyPolygon *)ctx->target;
	// query with rasterization
	if(ctx->use_grid){
		poly->rasterization(ctx->vpr);
		target->rasterization(ctx->vpr);
	}

	if(poly == target){
		return true;
	}
	// the minimum possible distance is larger than the threshold
	if(poly->getMBB()->distance(*target->getMBB(), ctx->geography)>ctx->within_distance){
        return true;
	}
	// the maximum possible distance is smaller than the threshold
	if(poly->getMBB()->max_distance(*target->getMBB(), ctx->geography)<=ctx->within_distance){
		ctx->found++;
        return true;
	}

	timeval start = get_cur_time();

	ctx->distance = poly->distance(target,ctx);
	double t1 = get_time_elapsed(start, true);
	ctx->distance = target->distance(poly,ctx);
	double t2 = get_time_elapsed(start, true);

	if(min(t1, t2)>0.1){

		bool judge[5];

		judge[0] = poly->area()>target->area();
		judge[1] = poly->area()/poly->getMBB()->area()<target->area()/target->getMBB()->area();
		judge[2] = poly->get_rastor()->get_step(false)>target->get_rastor()->get_step(false);
		judge[3] = poly->get_num_vertices()>target->get_num_vertices();
		judge[4] = tryluck(0.5);
		pthread_mutex_lock(&judge_t);
		total_checked++;
		for(int i=0;i<5;i++){
			judge_count[i] += ((t1<=t2)==judge[i]);
			if(judge[i]){
				right_total[i] += t1;
				wrong_total[i] += t2;
			}else{
				right_total[i] += t2;
				wrong_total[i] += t1;
			}
		}
		optimal += min(t1,t2);

		pthread_mutex_unlock(&judge_t);
	}



	ctx->found += ctx->distance <= ctx->within_distance;

	if(ctx->collect_latency){
		int nv = target->get_num_vertices();
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

	while(ctx->next_batch(1)){
		for(int i=ctx->index;i<ctx->index_end;i++){
			if(!tryluck(gctx->sample_rate)){
				continue;
			}
			struct timeval query_start = get_cur_time();
			ctx->target = (void *)(gctx->source_polygons[i]);
			box qb = gctx->source_polygons[i]->getMBB()->expand(gctx->within_distance, ctx->geography);
			tree.Search(qb.low, qb.high, MySearchCallback, (void *)ctx);
//			if(gctx->source_polygons[i]->getid()==8){
//				gctx->source_polygons[i]->print(false, false);
//			}
			ctx->object_checked.execution_time += get_time_elapsed(query_start);
			ctx->report_progress();
		}
	}
	ctx->merge_global();
	return NULL;
}


vector<MyPolygon *> load_binary_file(const char *path, query_context &ctx, bool sample){
	vector<MyPolygon *> polygons;
	if(!file_exist(path)){
		log("%s does not exist",path);
		return polygons;
	}
	struct timeval start = get_cur_time();
	ifstream infile;
	infile.open(path, ios::in | ios::binary);
	size_t off;
	size_t num_polygons = 0;
	infile.seekg(0, infile.end);
	//seek to the first polygon
	infile.seekg(-sizeof(size_t), infile.end);
	infile.read((char *)&num_polygons, sizeof(size_t));
	assert(num_polygons>0 && "the file should contain at least one polygon");
	size_t *offsets = new size_t[num_polygons];

	infile.seekg(-sizeof(size_t)*(num_polygons+1), infile.end);
	infile.read((char *)offsets, sizeof(size_t)*num_polygons);
	num_polygons = min(num_polygons, ctx.max_num_polygons);

	log("loading %ld polygon from %s",num_polygons,path);

	size_t num_edges = 0;
	size_t data_size = 0;
	for(size_t i=0;i<num_polygons;i++){
		if(i%100000==0){
			log("loaded %ld%%",i);
		}
		if(sample && !tryluck(ctx.sample_rate)){
			continue;
		}
		infile.seekg(offsets[i], infile.beg);
		MyPolygon *poly = MyPolygon::read_polygon_binary_file(infile);
		assert(poly);
		num_edges += poly->get_num_vertices();
		data_size += poly->get_data_size();
		poly->setid(i);
		polygons.push_back(MyPolygon::gen_box(*poly->getMBB()));
		delete poly;
	}
	infile.close();

	logt("loaded %ld polygons each with %ld edges %ld MB", start, polygons.size(),num_edges/polygons.size(),data_size/1024/1024);
	return polygons;
}


int main(int argc, char** argv) {



	query_context global_ctx;
	global_ctx = get_parameters(argc, argv);

	vector<MyPolygon *> boxes = load_binary_file(global_ctx.source_path.c_str(), global_ctx, false);;

	::dump_polygons_to_file(boxes, "all.boxes.dat");
	return 0;

	global_ctx.query_type = QueryType::within;
    global_ctx.geography = true;

	timeval start = get_cur_time();
//	global_ctx.max_num_polygons = 10000;
    global_ctx.source_polygons = MyPolygon::load_binary_file(global_ctx.source_path.c_str(), global_ctx);
	logt("loaded %d objects", start,global_ctx.source_polygons.size());

	preprocess(&global_ctx);
	logt("preprocess the source data", start);

	for(MyPolygon *p:global_ctx.source_polygons){
		tree.Insert(p->getMBB()->low, p->getMBB()->high, p);
	}
	logt("building R-Tree with %d nodes", start, global_ctx.source_polygons.size());

	// the target is also the source
	global_ctx.target_num = global_ctx.source_polygons.size();
	//global_ctx.target_num = 1;
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

	for(int i=0;i<5;i++){
		printf("%f\t%f\t%f\n", judge_count[i]*1.0/total_checked,
							   right_total[i]/optimal,
							   wrong_total[i]/optimal);
	}
	return 0;
}



