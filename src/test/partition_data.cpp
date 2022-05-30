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


class polygon_list{
public:
	int id = 0;
	vector<MyPolygon *> polygons;
	pthread_mutex_t lk;
	polygon_list(){
		pthread_mutex_init(&lk, NULL);
	}

};
// some shared parameters

RTree<polygon_list *, double, 2, double> tree;

bool MySearchCallback(polygon_list *polys, void* arg){
	query_context *ctx = (query_context *)arg;

	MyPolygon *target = (MyPolygon *)ctx->target;
	pthread_mutex_lock(&polys->lk);
    polys->polygons.push_back(target);
	pthread_mutex_unlock(&polys->lk);

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
			//gctx->source_polygons[i]->getMBB()->print();
			ctx->target = (void *)gctx->source_polygons[i];
			tree.Search(gctx->source_polygons[i]->getMBB()->low, gctx->source_polygons[i]->getMBB()->high, MySearchCallback, (void *)ctx);
			ctx->report_progress();
		}
	}
	ctx->merge_global();

	return NULL;
}



int main(int argc, char** argv) {

	MyMultiPolygon *mp = MyMultiPolygon::read_multipolygon();
	vector<MyPolygon *> polygons = mp->get_polygons();

	vector<polygon_list *> partitions;
	for(int i=0;i<polygons.size();i++){
		polygon_list *pos = new polygon_list();
		pos->id = i;
		partitions.push_back(pos);
		polygons[i]->setid(i);
		tree.Insert(polygons[i]->getMBB()->low, polygons[i]->getMBB()->high, partitions[i]);
	}

	query_context global_ctx;
	global_ctx = get_parameters(argc, argv);
	global_ctx.source_polygons = MyPolygon::load_binary_file(global_ctx.source_path.c_str(),global_ctx,true);
	global_ctx.target_num = global_ctx.source_polygons.size();
	timeval start = get_cur_time();

	// read all the points
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
	logt("total query",start);

	int total = 0;

	int max_one = 0;
	for(int i=0;i<partitions.size();i++){
		char path[256];
		sprintf(path, "part_polygons/%d.dat",i);

		dump_polygons_to_file(partitions[i]->polygons, path);

		log("%d %d",i,partitions[i]->polygons.size());
		if(partitions[i]->polygons.size()>partitions[max_one]->polygons.size()){
			max_one = i;
		}
	}
	log("%d %d",max_one,partitions[max_one]->polygons.size());

//	for(MyPolygon *p:partitions[max_one]->polygons){
//		p->print(false, false);
//	}
	//polygons[max_one]->print(false, false);

	return 0;
}



