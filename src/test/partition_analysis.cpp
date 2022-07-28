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


// functions for partitioning
bool assign(Tile *tile, void *arg){
	tile->insert((MyPolygon *)arg, false);
	return true;
}

void *partition_unit(void *arg){
	query_context *ctx = (query_context *)arg;
	vector<MyPolygon *> *objects = (vector<MyPolygon *> *)(ctx->target);
	RTree<Tile *, double, 2, double> *global_tree = (RTree<Tile *, double, 2, double> *)(ctx->target2);
	while(ctx->next_batch(100)){
		for(size_t i=ctx->index;i<ctx->index_end;i++){
			global_tree->Search((*objects)[i]->getMBB()->low, (*objects)[i]->getMBB()->high, assign, (void *)((*objects)[i]));
			ctx->report_progress();
		}
	}
	ctx->merge_global();
	return NULL;
}

void partition(vector<MyPolygon *> &objects, RTree<Tile *, double, 2, double> &global_tree, PARTITION_TYPE ptype){

	query_context global_ctx;
	pthread_t threads[global_ctx.num_threads];
	query_context ctx[global_ctx.num_threads];
	global_ctx.target_num = objects.size();
	global_ctx.report_prefix = "partitioning";

	for(int i=0;i<global_ctx.num_threads;i++){
		ctx[i] = global_ctx;
		ctx[i].thread_id = i;
		ctx[i].global_ctx = &global_ctx;
		ctx[i].target = (void *) &objects;
		ctx[i].target2 = (void *) &global_tree;
	}
	for(int i=0;i<global_ctx.num_threads;i++){
		pthread_create(&threads[i], NULL, partition_unit, (void *)&ctx[i]);
	}

	for(int i = 0; i < global_ctx.num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}

}

bool assign_target(Tile *tile, void *arg){
	query_context *ctx = (query_context *)arg;
	tile->insert_target((Point *)ctx->target);
	return true;
}

void *partition_target_unit(void *arg){
	query_context *ctx = (query_context *)arg;
	vector<Point *> *objects = (vector<Point *> *)(ctx->target);
	RTree<Tile *, double, 2, double> *global_tree = (RTree<Tile *, double, 2, double> *)(ctx->target2);
	while(ctx->next_batch(10)){
		for(size_t i=ctx->index;i<ctx->index_end;i++){
			Point *p = (*objects)[i];
			vector<Tile *> belongs;
			query_context tmpctx;
			tmpctx.target = (void *)p;
			global_tree->Search((double *)p, (double *)p, assign_target, (void *)&tmpctx);
			ctx->report_progress();
		}
	}
	ctx->merge_global();
	return NULL;
}

void partition_target(vector<Point *> &targets, RTree<Tile *, double, 2, double> &global_tree){

	query_context global_ctx;
	//global_ctx.num_threads = 1;
	pthread_t threads[global_ctx.num_threads];
	query_context ctx[global_ctx.num_threads];
	global_ctx.target_num = targets.size();
	global_ctx.report_prefix = "partitioning";

	for(int i=0;i<global_ctx.num_threads;i++){
		ctx[i] = global_ctx;
		ctx[i].thread_id = i;
		ctx[i].global_ctx = &global_ctx;
		ctx[i].target = (void *) &targets;
		ctx[i].target2 = (void *) &global_tree;
	}
	for(int i=0;i<global_ctx.num_threads;i++){
		pthread_create(&threads[i], NULL, partition_target_unit, (void *)&ctx[i]);
	}

	for(int i = 0; i < global_ctx.num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}

}

// functions for the query phase
void *local_unit(void *arg){
	query_context *ctx = (query_context *)arg;
	vector<Tile *> *tiles = (vector<Tile *> *)(ctx->target);
	for(size_t i=ctx->thread_id; i<ctx->target_num; i+= ctx->num_threads){
		Tile *tile = (*tiles)[i];
		tile->build_index();
		ctx->found += tile->conduct_query();
		ctx->next_batch(1);
		ctx->report_progress();
	}
	ctx->merge_global();
	return NULL;
}

size_t local(vector<Tile *> &tiles){

	query_context global_ctx;
	pthread_t threads[global_ctx.num_threads];
	query_context ctx[global_ctx.num_threads];
	global_ctx.target_num = tiles.size();
	global_ctx.report_prefix = "querying";

	for(int i=0;i<global_ctx.num_threads;i++){
		ctx[i] = global_ctx;
		ctx[i].thread_id = i;
		ctx[i].global_ctx = &global_ctx;
		ctx[i].target = (void *) &tiles;
	}

	for(int i=0;i<global_ctx.num_threads;i++){
		pthread_create(&threads[i], NULL, local_unit, (void *)&ctx[i]);
	}

	for(int i = 0; i < global_ctx.num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}
	return global_ctx.found;
}

typedef struct partition_stat_{
	double sample_rate = 1.0;
	size_t cardinality = 1;
	PARTITION_TYPE ptype;
	bool is_data_oriented = false;

	double stddev_r = 0; // for reference
	double stddev_t = 0; // for target
	double stddev_a = 0; // for all

	double boundary_rate = 0;
	double target_boundary_rate = 0.0;
	double found_rate = 0;
	size_t tile_num = 0;
	double genschema_time = 0.0;
	double partition_time = 0.0;
	double query_time = 0;
	double total_time = 0.0;

	void print(){
		printf("setup,%s, %f,%s,%ld,%ld,"
				"stddev,%f,%f,%f,"
				"stats,%f,%f,%f,"
				"time,%f,%f,%f,%f\n",
				is_data_oriented?"data":"space", sample_rate,partition_type_names[ptype], cardinality, tile_num,
				stddev_r, stddev_t, stddev_a,
				boundary_rate, target_boundary_rate, found_rate,
				genschema_time, partition_time, query_time, total_time);
		fflush(stdout);
	}
}partition_stat;

static partition_stat process(vector<vector<MyPolygon *>> &object_sets, vector<vector<Point *>> &target_sets, size_t card, PARTITION_TYPE ptype, bool data_oriented){

	assert(object_sets.size()==target_sets.size());
	partition_stat stat;
	stat.is_data_oriented = data_oriented;
	struct timeval start = get_cur_time();
	struct timeval entire_start = get_cur_time();

	// generate schema
	vector<Tile *> tiles;
	vector<RTree<Tile *, double, 2, double> *> global_trees;

	for(vector<MyPolygon *> &objects:object_sets){
		struct timeval cst = get_cur_time();
		vector<Tile *> tl = genschema(objects, card, ptype, data_oriented);
		logt("generating schema", cst);
		RTree<Tile *, double, 2, double> *global_tree = new RTree<Tile *, double, 2, double>();
		for(Tile *t:tl){
			global_tree->Insert(t->low, t->high, t);
		}
		logt("build global tree", cst);
		global_trees.push_back(global_tree);
		tiles.insert(tiles.end(), tl.begin(), tl.end());
		tl.clear();
	}
	stat.genschema_time = logt("%ld tiles are generated with %s partitioning algorithm",start, tiles.size(), partition_type_names[ptype]);

	// partitioning data
	if(!data_oriented && ptype!=HC){
		for(size_t i=0;i<object_sets.size();i++){
			partition(object_sets[i], *global_trees[i], ptype);
		}
	}
	for(size_t i=0;i<target_sets.size();i++){
		partition_target(target_sets[i], *global_trees[i]);
	}
	stat.partition_time = logt("partitioning data",start);

	// local processing
	size_t found = local(tiles);
	stat.query_time = logt("local processing get %ld pairs",start,found);
	stat.total_time = logt("entire process", entire_start);


	size_t total_num = 0;
	size_t total_target_num = 0;
	for(Tile *t:tiles){
		total_num += t->objects.size();
		total_target_num += t->targets.size();
	}

	stat.stddev_r = skewstdevratio(tiles, 0);
	stat.stddev_t = skewstdevratio(tiles, 1);
	stat.stddev_a = skewstdevratio(tiles, 2);

	size_t objnum = 0;
	size_t tgtnum = 0;
	for(vector<MyPolygon *> &objects:object_sets){
		objnum += objects.size();
	}
	for(vector<Point *> &targets:target_sets){
		tgtnum += targets.size();
	}
	//log("%ld\t%ld", total_num, objnum);
	stat.boundary_rate = 1.0*total_num/objnum;
	stat.target_boundary_rate = 1.0*total_target_num/tgtnum;
	stat.found_rate = 1.0*found/tgtnum;
	stat.tile_num = tiles.size();

//	size_t x_crossed = 0;
//	size_t y_crossed = 0;
//	for(Tile *tile:tiles){
//		for(MyPolygon *p:tile->objects){
//			x_crossed += (p->getMBB()->high[0]>tile->high[0]||p->getMBB()->low[0]<tile->low[0]);
//			y_crossed += (p->getMBB()->high[1]>tile->high[1]||p->getMBB()->low[1]<tile->low[1]);
//		}
//	}
//	log("%ld,%ld",x_crossed,y_crossed);

	// clear the partition schema for this round
	for(Tile *tile:tiles){
		delete tile;
	}
	tiles.clear();

	for(size_t i=0;i<target_sets.size();i++){
		delete global_trees[i];
	}

	return stat;
}

int main(int argc, char** argv) {

	string mbr_path;
	string point_path;

	size_t min_cardinality = 2000;
	size_t max_cardinality = 2000;
	bool fixed_cardinality = false;
	bool data_oriented = false;

	double min_sample_rate = 0.01;
	double max_sample_rate = 0.01;

	string ptype_str = "str";

	po::options_description desc("query usage");
	desc.add_options()
		("help,h", "produce help message")
		("fixed_cardinality,f", "fixed cardinality for all sampling rates")
		("data_oriented,d", "data-oriented partitioning")
		("source,s", po::value<string>(&mbr_path)->required(), "path to the source mbrs")
		("target,t", po::value<string>(&point_path), "path to the target points")

		("partition_type,p", po::value<string>(&ptype_str), "partition type should be one of: str|slc|qt|bsp|hc|fg, if not set, all the algorithms will be tested")

		("sample_rate,r", po::value<double>(&min_sample_rate), "the maximum sample rate (0.01 by default)")
		("max_sample_rate", po::value<double>(&max_sample_rate), "the maximum sample rate (0.01 by default)")

		("cardinality,c", po::value<size_t>(&min_cardinality), "the base number of objects per partition contain (2000 by default)")
		("max_cardinality", po::value<size_t>(&max_cardinality), "the max number of objects per partition contain (2000 by default)")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
		cout << desc << "\n";
		exit(0);
	}
	po::notify(vm);

	if(max_sample_rate<min_sample_rate||!vm.count("max_sample_rate")){
		max_sample_rate = min_sample_rate;
	}
	if(max_cardinality<min_cardinality || !vm.count("max_cardinality")){
		max_cardinality = min_cardinality;
	}
	fixed_cardinality = vm.count("fixed_cardinality");
	data_oriented = vm.count("data_oriented");

	PARTITION_TYPE start_type = STR;
	PARTITION_TYPE end_type = BSP;
	if(vm.count("partition_type")){
		PARTITION_TYPE ptype = ::parse_partition_type(ptype_str.c_str());
		start_type = ptype;
		end_type = ptype;
	}

	vector<vector<MyPolygon *>> object_sets;
	vector<vector<Point *>> target_sets;
	vector<Point *> raw_points;

	vector<string> objfiles;
	list_files(mbr_path.c_str(), objfiles);

	struct timeval start = get_cur_time();
	size_t numobj = 0;
	size_t numtgt = 0;
	if(vm.count("target")){
		assert(objfiles.size()==1);
		query_context ctx;
		ctx.sample_rate = max_sample_rate;
		vector<MyPolygon *> polygons = load_binary_file(objfiles[0].c_str(), ctx);
		numobj += polygons.size();
		object_sets.push_back(polygons);

		Point *points;
		size_t points_num = load_points_from_path(point_path.c_str(), &points);

		vector<Point *> point_vec;
		point_vec.resize(points_num);
		for(size_t i=0;i<points_num;i++){
			point_vec[i] = points+i;
		}
		vector<Point *> targets = sample<Point>(point_vec, max_sample_rate);
		raw_points.push_back(points);
		target_sets.push_back(targets);
		numtgt += targets.size();
	}else{
		for(string s:objfiles){
			query_context ctx;
			ctx.sample_rate = max_sample_rate;
			vector<MyPolygon *> polygons = load_binary_file(s.c_str(), ctx);
			numobj += polygons.size();
			object_sets.push_back(polygons);

			Point *points = new Point[polygons.size()];
			vector<Point *> targets;
			targets.resize(polygons.size());
			for(size_t i=0;i<polygons.size();i++){
				points[i] = polygons[i]->getMBB()->centroid();
				targets[i] = &points[i];
			}
			raw_points.push_back(points);
			target_sets.push_back(targets);
			numtgt += targets.size();
		}
	}

	logt("%ld objects and %ld targets are sampled from %ld files", start, numobj, numtgt, objfiles.size());
	// iterate the sampling rate
	for(double sample_rate = max_sample_rate;sample_rate*1.5>=min_sample_rate; sample_rate /= 2){
		vector<vector<MyPolygon *>> cur_sampled_objects;
		vector<vector<Point *>> cur_sampled_points;

		size_t cobjnum = 0;
		size_t ctgtnum = 0;
		for(size_t i=0;i<objfiles.size();i++){
			vector<MyPolygon *> cobj;
			vector<Point *> ctgt;
			if(sample_rate == max_sample_rate){
				cobj.insert(cobj.begin(), object_sets[i].begin(), object_sets[i].end());
				ctgt.insert(ctgt.begin(), target_sets[i].begin(), target_sets[i].end());
			}else{
				cobj = sample<MyPolygon>(object_sets[i], sample_rate/max_sample_rate);
				ctgt = sample<Point>(target_sets[i], sample_rate/max_sample_rate);
			}
			cobjnum += cobj.size();
			ctgtnum += ctgt.size();
			cur_sampled_objects.push_back(cobj);
			cur_sampled_points.push_back(ctgt);
		}

		logt("resampled %ld MBRs and %ld points with sample rate %f", start,cobjnum,ctgtnum, sample_rate);

		for(int pt=(int)start_type;pt<=(int)end_type;pt++){
			vector<double> num_tiles;
			//vector<double> stddevs;
			vector<double> boundary;
			vector<double> found;
			// varying the cardinality
			for(size_t card=min_cardinality; (card/1.5)<=max_cardinality;card*=2){
				size_t real_card = std::max((int)(card*sample_rate), 1);
				PARTITION_TYPE ptype = (PARTITION_TYPE)pt;
				log("processing %f\t%s\t%ld",sample_rate,partition_type_names[ptype],real_card);
				partition_stat stat = process(cur_sampled_objects, cur_sampled_points, real_card, ptype, data_oriented);
				stat.sample_rate = sample_rate,
				stat.ptype = ptype;
				stat.cardinality = real_card;
				stat.print();
				logt("complete",start);
				num_tiles.push_back(stat.tile_num);
				//stddevs.push_back(stat.stddev);
				boundary.push_back(stat.boundary_rate);
				found.push_back(stat.found_rate);
			}
			if(num_tiles.size()>1){
				double a,b;
				//linear_regression(num_tiles,stddevs,a,b);
				//printf("%f,%.10f,%.10f",sample_rate,a,b);
				linear_regression(num_tiles,boundary,a,b);
				printf(",%.10f,%.10f",a,b);
				linear_regression(num_tiles,found,a,b);
				printf(",%.10f,%.10f\n",a,b);
			}
		}
		for(vector<MyPolygon *> &ps:cur_sampled_objects){
			ps.clear();
		}
		for(vector<Point *> &ps:cur_sampled_points){
			ps.clear();
		}
		cur_sampled_objects.clear();
		cur_sampled_points.clear();
	}

	for(vector<MyPolygon *> &ps:object_sets){
#pragma omp parallel for num_threads(2*get_num_threads()-1)
		for(MyPolygon *p:ps){
			delete p;
		}
		ps.clear();
	}
	object_sets.clear();
	for(vector<Point *> &ps:target_sets){
		ps.clear();
	}

	for(Point *p:raw_points){
		delete []p;
	}
	raw_points.clear();
	return 0;
}
