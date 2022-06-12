/*
 * partition_analysis.cpp
 *
 *  Created on: Jun 8, 2022
 *      Author: teng
 */

#include "partition.h"

namespace po = boost::program_options;

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

	string mbr_path;
	string point_path;

	size_t min_cardinality = 2000;
	size_t max_cardinality = 2000;

	double min_sample_rate = 0.01;
	double max_sample_rate = 0.01;

	int sample_rounds = 1;

	string ptype_str = "str";

	po::options_description desc("query usage");
	desc.add_options()
		("help,h", "produce help message")
		("source,s", po::value<string>(&mbr_path)->required(), "path to the source mbrs")
		("target,t", po::value<string>(&point_path)->required(), "path to the target points")

		("partition_type,p", po::value<string>(&ptype_str), "partition type should be one of: str|slc|bos|qt|bsp|hc|fg, if not set, all the algorithms will be tested")

		("sample_rate,r", po::value<double>(&min_sample_rate), "the maximum sample rate (0.01 by default)")
		("max_sample_rate", po::value<double>(&max_sample_rate), "the maximum sample rate (0.01 by default)")

		("cardinality,c", po::value<size_t>(&min_cardinality), "the base number of objects per partition contain (2000 by default)")
		("max_cardinality", po::value<size_t>(&max_cardinality), "the max number of objects per partition contain (2000 by default)")

		("sample_rounds", po::value<int>(&sample_rounds), "the number of rounds the tests should be repeated")

		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
		cout << desc << "\n";
		exit(0);
	}
	po::notify(vm);

	if(max_sample_rate<min_sample_rate){
		max_sample_rate = min_sample_rate;
	}
	if(max_cardinality<min_cardinality){
		max_cardinality = min_cardinality;
	}

	PARTITION_TYPE start_type = STR;
	PARTITION_TYPE end_type = BSP;
	if(vm.count("partition_type")){
		PARTITION_TYPE ptype = ::parse_partition_type(ptype_str.c_str());
		start_type = ptype;
		end_type = ptype;
	}

	struct timeval start = get_cur_time();
	box *boxes;
	size_t box_num = load_boxes_from_file(mbr_path.c_str(), &boxes);
	vector<box> box_vec(boxes, boxes+box_num);
	logt("%ld objects are loaded",start, box_num);

	Point *points;
	size_t points_num = load_points_from_path(point_path.c_str(), &points);
	vector<Point> point_vec(points, points+points_num);

	logt("%ld points are loaded", start, points_num);

	for(int sr=0;sr<sample_rounds;sr++){
		// iterate the sampling rate
		for(double sample_rate = min_sample_rate;sample_rate/2.0<=max_sample_rate; sample_rate *= 2){
			// sampling data
			vector<box *> objects = sample<box>(box_vec, sample_rate);
			logt("%ld objects are sampled with sample rate %f",start, objects.size(),sample_rate);

			vector<Point *> targets = sample<Point>(point_vec, sample_rate);
			logt("%ld targets are sampled with sample rate %f",start, targets.size(),sample_rate);

			for(int pt=(int)start_type;pt<=(int)end_type;pt++){
				PARTITION_TYPE ptype = (PARTITION_TYPE)pt;

				// iterate the partitioning number
				for(size_t card=min_cardinality; (card/2)<=max_cardinality;card*=2){
					// generate schema
					vector<Tile *> tiles = genschema(objects, card, ptype);

					RTree<Tile *, double, 2, double> tree;
					for(Tile *t:tiles){
						tree.Insert(t->low, t->high, t);
					}
					logt("%ld tiles are generated with %s partitioning algorithm",start, tiles.size(), partition_type_names[ptype]);

					// partitioning data
					partition(objects, tree);
					size_t total_num = 0;
					for(Tile *tile:tiles){
						total_num += tile->get_objnum();
					}
					double partition_time = get_time_elapsed(start);
					logt("partitioning data",start);

					// conduct query
					size_t found = query(targets, tree);
					double query_time = get_time_elapsed(start);
					logt("querying data",start);
					printf("%s,%d,%f,%ld,%ld,%ld,%ld,%ld,%f,%f\n",
							partition_type_names[ptype],sr,
							sample_rate, card,
							objects.size(), total_num,
							targets.size(), found,
							partition_time, query_time);
					// clear the partition schema for this round
					for(Tile *tile:tiles){
						delete tile;
					}
					tiles.clear();
				}
			}
			objects.clear();
			targets.clear();
		}
	}

	delete []boxes;
	delete []points;
	return 0;
}
