/*
 * partition_analysis.cpp
 *
 *  Created on: Jun 8, 2022
 *      Author: teng
 */

#include "partition.h"

// functions for sampling

template <class O>
void *sample_unit(void *arg){
	query_context *ctx = (query_context *)arg;
	vector<O> *original = (vector<O> *)(ctx->target);
	vector<O *> *result = (vector<O *> *)(ctx->target2);
	while(ctx->next_batch(1000)){
		for(size_t i=ctx->index;i<ctx->index_end;i++){
			if(tryluck(ctx->global_ctx->sample_rate)){
				ctx->global_ctx->lock();
				result->push_back(&(*original)[i]);
				ctx->global_ctx->unlock();
			}
			ctx->report_progress();
		}
	}
	ctx->merge_global();
	return NULL;
}

template <class O>
vector<O *> sample(vector<O> &original, double sample_rate){

	vector<O *> result;
	query_context global_ctx;
	pthread_t threads[global_ctx.num_threads];
	query_context ctx[global_ctx.num_threads];
	global_ctx.target_num = original.size();
	global_ctx.sample_rate = sample_rate;
	for(int i=0;i<global_ctx.num_threads;i++){
		ctx[i] = global_ctx;
		ctx[i].thread_id = i;
		ctx[i].global_ctx = &global_ctx;
		ctx[i].target = (void *) &original;
		ctx[i].target2 = (void *) &result;
	}
	for(int i=0;i<global_ctx.num_threads;i++){
		pthread_create(&threads[i], NULL, sample_unit<O>, (void *)&ctx[i]);
	}

	for(int i = 0; i < global_ctx.num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}

	return result;
}

// functions for partitioning
bool assign(Tile *tile, void *arg){
	tile->insert((box *)arg, (box *)arg);
	return true;
}

void *partition_unit(void *arg){
	query_context *ctx = (query_context *)arg;
	vector<box *> *objects = (vector<box *> *)(ctx->target);
	RTree<Tile *, double, 2, double> *tree= (RTree<Tile *, double, 2, double> *)(ctx->target2);
	while(ctx->next_batch(100)){
		for(size_t i=ctx->index;i<ctx->index_end;i++){
			tree->Search((*objects)[i]->low, (*objects)[i]->high, assign, (void *)((*objects)[i]));
			ctx->report_progress();
		}
	}
	ctx->merge_global();
	return NULL;
}

void partition(vector<box *> &objects, RTree<Tile *, double, 2, double> &tree){

	query_context global_ctx;
	pthread_t threads[global_ctx.num_threads];
	query_context ctx[global_ctx.num_threads];
	global_ctx.target_num = objects.size();
	for(int i=0;i<global_ctx.num_threads;i++){
		ctx[i] = global_ctx;
		ctx[i].thread_id = i;
		ctx[i].global_ctx = &global_ctx;
		ctx[i].target = (void *) &objects;
		ctx[i].target2 = (void *) &tree;
	}
	for(int i=0;i<global_ctx.num_threads;i++){
		pthread_create(&threads[i], NULL, partition_unit, (void *)&ctx[i]);
	}

	for(int i = 0; i < global_ctx.num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}

}

// function for the query phase
bool search(Tile *tile, void *arg){
	query_context *ctx = (query_context *)arg;
	ctx->found += tile->lookup_count((Point *)ctx->target);
	return true;
}

void *query_unit(void *arg){
	query_context *ctx = (query_context *)arg;
	vector<Point *> *objects = (vector<Point *> *)(ctx->target);
	RTree<Tile *, double, 2, double> *tree= (RTree<Tile *, double, 2, double> *)(ctx->target2);
	while(ctx->next_batch(100)){
		for(size_t i=ctx->index;i<ctx->index_end;i++){
			ctx->target = (void *) (*objects)[i];
			tree->Search((double *)(*objects)[i], (double *)(*objects)[i], search, (void *)ctx);
			ctx->report_progress();
		}
	}
	ctx->merge_global();
	return NULL;
}

size_t query(vector<Point *> &objects, RTree<Tile *, double, 2, double> &tree){

	query_context global_ctx;
	pthread_t threads[global_ctx.num_threads];
	query_context ctx[global_ctx.num_threads];
	global_ctx.target_num = objects.size();
	for(int i=0;i<global_ctx.num_threads;i++){
		ctx[i] = global_ctx;
		ctx[i].thread_id = i;
		ctx[i].global_ctx = &global_ctx;
		ctx[i].target = (void *) &objects;
		ctx[i].target2 = (void *) &tree;
	}
	for(int i=0;i<global_ctx.num_threads;i++){
		pthread_create(&threads[i], NULL, query_unit, (void *)&ctx[i]);
	}

	for(int i = 0; i < global_ctx.num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}
	return global_ctx.found;

}

int main(int argc, char** argv) {

	const char *mbr_path = argv[1];
	const char *point_path = argv[2];
	const size_t partition_num = atoi(argv[3]);
	const PARTITION_TYPE ptype = parse_partition_type(argv[4]);
	const double max_sample_rate = 0.01;
	const size_t min_sample_ratio = 1;
	const size_t max_partition_ratio = 1;

	struct timeval start = get_cur_time();
	box *boxes;
	size_t box_num = load_boxes_from_file(mbr_path, &boxes);
	vector<box> box_vec(boxes, boxes+box_num);
	logt("%ld objects are loaded",start, box_num);

	Point *points;
	size_t points_num = load_points_from_path(point_path, &points);
	vector<Point> point_vec(points, points+points_num);

	logt("%ld points are loaded", start, points_num);
	// iterate the sampling rate
	for(size_t s=1;s<=min_sample_ratio;s*=2){
		// sampling data
		double sample_rate = max_sample_rate*s/min_sample_ratio;
		vector<box *> objects = sample<box>(box_vec, sample_rate);
		logt("%ld objects are sampled",start, objects.size());

		vector<Point *> targets = sample<Point>(point_vec, sample_rate);
		logt("%ld targets are sampled",start, targets.size());

		// iterate the partitioning number
		for(size_t pr=1;pr<=max_partition_ratio;pr*=2){
			// generate schema
			size_t true_partition_num = partition_num*pr;
			vector<Tile *> tiles = genschema(objects, true_partition_num, ptype);
			RTree<Tile *, double, 2, double> tree;
			for(Tile *t:tiles){
				tree.Insert(t->low, t->high, t);
			}
			logt("%ld tiles are generated",start, tiles.size());

			// partitioning data
			partition(objects, tree);
			size_t total_num = 0;
			for(Tile *tile:tiles){
				total_num += tile->get_objnum();
			}
			logt("partitioning data",start);

			// conduct query
			size_t found = query(targets, tree);
			logt("%.4f boundary %.4f (%ld/%ld) average hit", start, 100.0*(total_num-objects.size())/objects.size(), 1.0*found/targets.size(),found, targets.size());

			// clear the partition schema for this round
			for(Tile *tile:tiles){
				delete tile;
			}
			tiles.clear();
		}
		objects.clear();
		targets.clear();
	}


	delete []boxes;
	delete []points;
	return 0;
}
