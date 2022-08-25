/*
 * partition_analysis.cpp
 *
 *  Created on: Jun 8, 2022
 *      Author: teng
 */

#include "partition.h"
#include "query_context.h"

namespace po = boost::program_options;

// functions for sampling
box china(73.47431,3.393568,134.86609,53.65586);

template <class O>
void *sample_unit(void *arg){
	query_context *ctx = (query_context *)arg;
	vector<O *> *original = (vector<O *> *)(ctx->target);
	vector<O *> *result = (vector<O *> *)(ctx->target2);
	while(ctx->next_batch(1000)){
		for(size_t i=ctx->index;i<ctx->index_end;i++){
			if(tryluck(ctx->global_ctx->sample_rate)){
				ctx->global_ctx->lock();
				result->push_back((*original)[i]);
				ctx->global_ctx->unlock();
			}
			ctx->report_progress();
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
	global_ctx.report_prefix = "sample";
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

	string data_path;
	size_t cardinality = 100;
	double sample_rate = 1.0;
	string ptype_str = "str";

	po::options_description desc("query usage");
	desc.add_options()
		("help,h", "produce help message")
		("source,s", po::value<string>(&data_path)->required(), "path to the source")
		("sample_rate,r", po::value<double>(&sample_rate), "sample rate")
		("cardinality,c", po::value<size_t>(&cardinality), "the number of objects per partition contain")
		("partition_type,t", po::value<string>(&ptype_str), "partition type should be one of: str|slc|qt|bsp|hc|fg")
		("print,p", "print the generated schema")
		("is_points", "the input is points")
		("data_oriented,d", "data oriented partitioning")
		("bottom_up", "generate the schema with bottom up method")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
		cout << desc << "\n";
		exit(0);
	}
	po::notify(vm);
	struct timeval start = get_cur_time();
	vector<MyPolygon *> all_objects;
	if(vm.count("is_points")){
		Point *points = NULL;
		size_t box_num = load_points_from_path(data_path.c_str(), &points);
		for(size_t i=0;i<box_num;i++){
			box b;
			b.low[0] = points[i].x;
			b.low[1] = points[i].y;
			b.high[0] = points[i].x;
			b.high[1] = points[i].y;
			all_objects.push_back(MyPolygon::gen_box(b));
		}
		delete []points;
	}else{
		query_context ctx;
		all_objects = load_binary_file(data_path.c_str(), ctx);
	}

	logt("%ld objects are loaded",start, all_objects.size());
	vector<MyPolygon *> objects;
	// sampling data
	if(sample_rate >= 1.0){
		objects.insert(objects.end(), all_objects.begin(), all_objects.end());
	}else{
		objects = sample<MyPolygon>(all_objects, sample_rate);
		logt("%ld objects are sampled",start, objects.size());
	}

	PARTITION_TYPE start_type = STR;
	PARTITION_TYPE end_type = BSP;

	if(vm.count("partition_type")){
		PARTITION_TYPE ptype = ::parse_partition_type(ptype_str.c_str());
		start_type = ptype;
		end_type = ptype;
	}
	for(int pt=start_type;pt<=end_type;pt++){
		vector<Tile *> tiles;
		if(!vm.count("bottom_up")){
			tiles = genschema(objects, std::max(1, (int)(sample_rate*cardinality)), (PARTITION_TYPE)pt, vm.count("data_oriented"));
		}else{
			tiles = genschema_st(objects, std::max(1, (int)(sample_rate*cardinality)), (PARTITION_TYPE)pt, vm.count("data_oriented"));
		}
		if(vm.count("print")){
			print_tiles(tiles);
		}
		logt("%ld tiles are generated with %s partitioning algorithm",start, tiles.size(), partition_type_names[pt]);
		// clear the partition schema for this round
		for(Tile *tile:tiles){
			delete tile;
		}
		tiles.clear();
	}

	for(MyPolygon *p:all_objects){
		delete p;
	}
	objects.clear();
	return 0;
}
