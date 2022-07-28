/*
 * partition_stats.cpp
 *
 *  Created on: Jul 25, 2022
 *      Author: teng
 */



#include "partition.h"

namespace po = boost::program_options;
// functions for sampling
template <class O>
void *sample_unit(void *arg){
	query_context *ctx = (query_context *)arg;
	vector<O *> *original = (vector<O *> *)(ctx->target);
	vector<O *> *result = (vector<O *> *)(ctx->target2);
	vector<O *> tmp;
	while(ctx->next_batch(10000)){
		for(size_t i=ctx->index;i<ctx->index_end;i++){
			if(tryluck(ctx->global_ctx->sample_rate)){
				tmp.push_back((*original)[i]);
			}
			ctx->report_progress();
		}
		if(tmp.size()>0){
			ctx->global_ctx->lock();
			result->insert(result->end(), tmp.begin(), tmp.end());
			ctx->global_ctx->unlock();
			tmp.clear();
		}
	}
	ctx->merge_global();
	return NULL;
}

template <class O>
vector<O *> sample(vector<O *> &original, double sample_rate){

	vector<O *> result;
	query_context global_ctx;
	pthread_t threads[global_ctx.num_threads];
	query_context ctx[global_ctx.num_threads];
	global_ctx.target_num = original.size();
	global_ctx.sample_rate = sample_rate;
	global_ctx.report_prefix = "sampling";
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

int main(int argc, char** argv) {

	string mbr_path;
	string point_path;

	size_t max_num = LONG_MAX;

	po::options_description desc("query usage");
	desc.add_options()
		("help,h", "produce help message")
		//("source,s", po::value<string>(&mbr_path)->required(), "path to the source mbrs")
		("target,t", po::value<string>(&point_path), "path to the target points")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
		cout << desc << "\n";
		exit(0);
	}
	po::notify(vm);

	struct timeval start = get_cur_time();

	Point *points;
	size_t points_num = load_points_from_path(point_path.c_str(), &points);

	vector<Point *> point_vec;
	point_vec.resize(points_num);
	for(size_t i=0;i<points_num;i++){
		point_vec[i] = points+i;
	}
	logt("loading points", start);

	for(double sr = 0.1;sr>0.0001;sr/=2){
		start = get_cur_time();
		vector<Point *> targets = sample<Point>(point_vec, sr);
		logt("sampled %ld objects", start, targets.size());
		RTree<Point *, double, 2, double> tree;
		start = get_cur_time();
		for(Point *p:targets){
			tree.Insert((double *)p, (double *)p, p);
		}

		printf("%f\t%f\n",sr, get_time_elapsed(start,true));
		targets.clear();

	}

	return 0;
}

