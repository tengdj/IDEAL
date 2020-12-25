/*
 * query_context.cpp
 *
 *  Created on: Sep 3, 2020
 *      Author: teng
 */
#include "MyPolygon.h"

query_context::query_context(){
	num_threads = get_num_threads();
	pthread_mutex_init(&lock, NULL);
}
query_context::~query_context(){
	candidates.clear();
//	for(MyPolygon *p:source_polygons){
//		delete p;
//	}
//	for(MyPolygon *p:target_polygons){
//		delete p;
//	}
//	if(points){
//		delete []points;
//	}

}

query_context::query_context(query_context &t){
	thread_id = t.thread_id;
	num_threads = t.num_threads;
	vpr = t.vpr;
	vpr_end = t.vpr_end;
	use_grid = t.use_grid;
	use_qtree = t.use_qtree;
	use_mer = t.use_mer;
	mer_sample_round = t.mer_sample_round;
	use_convex_hull = t.use_convex_hull;
	perform_refine = t.perform_refine;
	gpu = t.gpu;
	sample_rate = t.sample_rate;
	small_threshold = t.small_threshold;
	big_threshold = t.big_threshold;
	sort_polygons = t.sort_polygons;
	report_gap = t.report_gap;
	distance_buffer_size = t.distance_buffer_size;
	source_path = t.source_path;
	target_path = t.target_path;
	query_type = t.query_type;
	collect_latency = t.collect_latency;
	pthread_mutex_init(&lock, NULL);
}

query_context& query_context::operator=(query_context const &t){
	thread_id = t.thread_id;
	num_threads = t.num_threads;
	vpr = t.vpr;
	vpr_end = t.vpr_end;
	use_grid = t.use_grid;
	use_qtree = t.use_qtree;
	use_mer = t.use_mer;
	mer_sample_round = t.mer_sample_round;
	use_convex_hull = t.use_convex_hull;
	perform_refine = t.perform_refine;
	gpu = t.gpu;
	sample_rate = t.sample_rate;
	small_threshold = t.small_threshold;
	big_threshold = t.big_threshold;
	sort_polygons = t.sort_polygons;
	report_gap = t.report_gap;
	distance_buffer_size = t.distance_buffer_size;
	source_path = t.source_path;
	target_path = t.target_path;
    query_type = t.query_type;
    collect_latency = t.collect_latency;
    pthread_mutex_init(&lock, NULL);
	return *this;

}

query_context query_context::operator+ (query_context const &obj) {
	query_context res = *this;

	return res;
}
void query_context::report_latency(int num_v, double lt){
	if(vertex_number.find(num_v)==vertex_number.end()){
		vertex_number[num_v] = 1;
		latency[num_v] = lt;
	}else{
		vertex_number[num_v] = vertex_number[num_v]+1;
		latency[num_v] = latency[num_v]+lt;
	}
}

void query_context::load_points(){
	struct timeval start = get_cur_time();
	long fsize = file_size(target_path.c_str());
	if(fsize<=0){
		log("%s is empty",target_path.c_str());
		exit(0);
	}else{
		log("size of %s is %ld",target_path.c_str(),fsize);
	}
	target_num = fsize/(2*sizeof(double));

	points = new double[target_num*2];
	ifstream infile(target_path.c_str(), ios::in | ios::binary);
	infile.read((char *)points, fsize);
	infile.close();
	logt("loaded %ld points", start,target_num);
}

void query_context::report_progress(){
	if(++query_count==1000){
		pthread_mutex_lock(&global_ctx->lock);
		global_ctx->query_count += query_count;
		if(global_ctx->query_count%global_ctx->report_gap==0){
			log("processed %d (%.2f\%)",global_ctx->query_count,(double)global_ctx->query_count*100/global_ctx->target_num);
		}
		query_count = 0;
		pthread_mutex_unlock(&global_ctx->lock);
	}
}

void query_context::merge_global(){
	pthread_mutex_lock(&global_ctx->lock);
	global_ctx->found += found;
	global_ctx->query_count += query_count;
	global_ctx->checked_count += checked_count;
	global_ctx->pixel_checked += pixel_checked;
	global_ctx->border_checked += border_checked;
	global_ctx->edges_checked += edges_checked;
	global_ctx->pixel_check_time += pixel_check_time;
	global_ctx->edges_check_time += edges_check_time;
	global_ctx->check_time += check_time;

	for(auto &it :vertex_number){
		const double lt = latency.at(it.first);
		if(global_ctx->vertex_number.find(it.first)!=global_ctx->vertex_number.end()){
			global_ctx->vertex_number[it.first] = global_ctx->vertex_number[it.first]+it.second;
			global_ctx->latency[it.first] = global_ctx->latency[it.first]+lt;
		}else{
			global_ctx->vertex_number[it.first] = it.second;
			global_ctx->latency[it.first] = lt;
		}
	}
	pthread_mutex_unlock(&global_ctx->lock);
}

bool query_context::next_batch(int batch_num){
	pthread_mutex_lock(&global_ctx->lock);
	if(global_ctx->index==global_ctx->target_num){
		pthread_mutex_unlock(&global_ctx->lock);
		return false;
	}
	index = global_ctx->index;
	if(index+batch_num>global_ctx->target_num){
		index_end = global_ctx->target_num;
	}else {
		index_end = index+batch_num;
	}
	global_ctx->index = index_end;
	pthread_mutex_unlock(&global_ctx->lock);
	return true;
}

void query_context::print_stats(){

	log("query count:\t%ld",query_count);
	log("checked count:\t%ld",checked_count);
	log("found count:\t%ld",found);

	if(checked_count>0){
		log("pixel/checked:\t%f",(double)pixel_checked/checked_count);
		log("border/checked:\t%f",(double)border_checked/checked_count);
		log("edges/checked:\t%f",(double)edges_checked/checked_count);
	}
	if(border_checked>0){
		log("edges/border:\t%f",(double)edges_checked/border_checked);
	}
	if(pixel_check_time>0){
		log("latency/pixel:\t%f",pixel_check_time/pixel_checked);
	}
	if(edges_check_time>0){
		log("latency/edge:\t%f",edges_check_time/edges_checked);
	}
	if(check_time>0){
		log("latency/other:\t%f",(check_time-edges_check_time-pixel_check_time)/checked_count);
	}

	if(collect_latency){
		for(auto it:vertex_number){
			cout<<it.first<<"\t"<<latency[it.first]/it.second<<endl;
		}
	}
}


query_context get_parameters(int argc, char **argv){
	query_context global_ctx;

	po::options_description desc("query usage");
	desc.add_options()
		("help,h", "produce help message")
		("rasterize,r", "partition with rasterization")
		("qtree,q", "partition with qtree")
		("raster_only", "query with raster only")
		("convex_hull", "use convex hall for filtering")
		("mer", "use maximum enclosed rectangle")
		("mer_sample_round", "how many rounds of sampling needed for MER generating")

		("source,s", po::value<string>(&global_ctx.source_path), "path to the source")
		("target,t", po::value<string>(&global_ctx.target_path), "path to the target")
		("threads,n", po::value<int>(&global_ctx.num_threads), "number of threads")
		("vpr,v", po::value<int>(&global_ctx.vpr), "number of vertices per raster")
		("vpr_end", po::value<int>(&global_ctx.vpr_end), "number of vertices per raster")
		("big_threshold,b", po::value<int>(&global_ctx.big_threshold), "up threshold for complex polygon")
		("small_threshold", po::value<int>(&global_ctx.small_threshold), "low threshold for complex polygon")
		("sample_rate", po::value<float>(&global_ctx.sample_rate), "sample rate")
		("latency,l","collect the latency information")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
		cout << desc << "\n";
		exit(0);
	}
	po::notify(vm);

	if(!vm.count("source")||!vm.count("target")){
		cout << desc << "\n";
		exit(0);
	}
	global_ctx.use_grid = vm.count("rasterize");
	global_ctx.use_qtree = vm.count("qtree");
	global_ctx.perform_refine = !vm.count("raster_only");
	global_ctx.collect_latency = vm.count("latency");
	global_ctx.use_convex_hull = vm.count("convex_hull");
	global_ctx.use_mer = vm.count("mer");


	assert(!(global_ctx.use_grid&&global_ctx.use_qtree));
	return global_ctx;
}
