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
		("partition_type,t", po::value<string>(&ptype_str), "partition type should be one of: str|slc|bos|qt|bsp|hc|fg")
		("print,p", "print the generated schema")

		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
		cout << desc << "\n";
		exit(0);
	}
	po::notify(vm);
	struct timeval start = get_cur_time();
	box *boxes;
	size_t box_num = load_boxes_from_file(data_path.c_str(), &boxes);
	logt("%ld objects are loaded",start, box_num);

	// sampling data
	vector<box> box_vec(boxes, boxes+box_num);
	vector<box *> objects = sample<box>(box_vec, sample_rate);
	logt("%ld objects are sampled",start, objects.size());

	PARTITION_TYPE start_type = STR;
	PARTITION_TYPE end_type = BSP;

	if(vm.count("partition_type")){
		PARTITION_TYPE ptype = ::parse_partition_type(ptype_str.c_str());
		start_type = ptype;
		end_type = ptype;
	}
	for(int pt=start_type;pt<=end_type;pt++){
		vector<Tile *> tiles = genschema(objects, std::max(1, (int)(sample_rate*cardinality)), (PARTITION_TYPE)pt);
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

	box_vec.clear();
	objects.clear();
	delete []boxes;
	return 0;
}
