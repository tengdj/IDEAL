///*
// * partstat.cpp
// *
// *  Created on: Jun 4, 2020
// *      Author: teng
// */
///*
// * Parser.cpp
// *
// *  Created on: May 9, 2020
// *      Author: teng
// */
//
//#include "../geometry/MyPolygon.h"
//#include <fstream>
//#include "../index/RTree.h"
//#include <queue>
//#include <boost/program_options.hpp>
//#include "cuda/mygpu.h"
//
//namespace po = boost::program_options;
//using namespace std;
//
//int element_size = 10;
//
//// some shared parameters
//pthread_mutex_t poly_lock;
//bool stop = false;
//
//pthread_mutex_t report_lock;
//long query_count = 0;
//
//MyPolygon *global_space = NULL;
//
//
//queue<queue_element *> shared_queue;
//int num_queries = 10;
//
//class global_info{
//public:
//	int borders = 0;
//	int ins = 0;
//	int outs = 0;
//	int crosses_count = 0;
//	int two_crosses_count = 0;
//	int multiple_crosses_count = 0;
//	int line_count = 0;
//	size_t num_pixels = 0;
//	double partition_time = 0;
//	double query_time = 0;
//	double nopartition_query_time = 0;
//	double distance_query_time = 0;
//	double nopartition_distance_query_time = 0;
//	void merge(global_info &target){
//		borders += target.borders;
//		ins += target.ins;
//		outs += target.outs;
//		crosses_count += target.crosses_count;
//		line_count += target.line_count;
//		num_pixels += target.num_pixels;
//		partition_time += target.partition_time;
//		query_time += target.query_time/num_queries;
//		nopartition_query_time += target.nopartition_query_time/num_queries;
//
//		distance_query_time += target.distance_query_time/num_queries;
//		nopartition_distance_query_time += target.nopartition_distance_query_time/num_queries;
//		this->two_crosses_count += target.two_crosses_count;
//		this->multiple_crosses_count += target.multiple_crosses_count;
//	}
//};
//
//class partition_context{
//public:
//	int vpr = 10;
//	int query_iteration = 0;
//	vector<int> vprs;
//	vector<double> flatneses;
//	vector<global_info> infos;
//	size_t total_num_vertices = 0;
//};
//
//
//
//
//void *partition(void *args){
//	partition_context *ctx = (partition_context *)args;
//	vector<queue_element *> elems;
//	int totalcount = ctx->vprs.size();
//	global_info local[totalcount];
//
//	size_t num_vertices = 0;
//	query_context qctx;
//	while(!stop||shared_queue.size()>0){
//		queue_element *elem = NULL;
//		pthread_mutex_lock(&poly_lock);
//		if(!shared_queue.empty()){
//			elem = shared_queue.front();
//			shared_queue.pop();
//		}
//		pthread_mutex_unlock(&poly_lock);
//		// queue is empty but not stopped
//		if(!elem){
//			usleep(20);
//			continue;
//		}
//		for(MyPolygon *poly:elem->polys){
//			num_vertices += poly->get_num_vertices();
//			vector<Point> targets = poly->generate_test_points(num_queries);
//			vector<Point> distance_targets = global_space->generate_test_points(num_queries);
//
//			for(int k=0;k<totalcount;k++){
//				timeval start = get_cur_time();
//				vector<vector<Pixel>> partitions = poly->partition(ctx->vprs[k]);
//				local[k].partition_time += get_time_elapsed(start,true);
//				for(int t=0;t<num_queries;t++){
//					qctx.use_grid = true;
//					poly->contain(targets[t],&qctx);
//					local[k].query_time += get_time_elapsed(start,true);
//					qctx.use_grid = false;
//					poly->contain(targets[t],&qctx);
//					local[k].nopartition_query_time += get_time_elapsed(start,true);
//					qctx.use_grid = true;
//					poly->distance(distance_targets[t],&qctx);
//					local[k].distance_query_time += get_time_elapsed(start,true);
//					qctx.use_grid = false;
//					poly->distance(distance_targets[t],&qctx);
//					local[k].nopartition_distance_query_time += get_time_elapsed(start,true);
//				}
//				local[k].line_count += partitions.size();
//				local[k].line_count += partitions[0].size();
//				for(int i=0;i<partitions.size();i++){
//					for(int j=0;j<partitions[0].size();j++){
//						if(partitions[i][j].status==BORDER){
//							local[k].borders++;
//						}else if(partitions[i][j].status==IN){
//							local[k].ins++;
//						}else if(partitions[i][j].status==OUT){
//							local[k].outs++;
//						}
//						local[k].crosses_count+=partitions[i][j].crosses.size();
//						local[k].num_pixels++;
//						if(partitions[i][j].crosses.size()==2){
//							local[k].two_crosses_count++;
//						}else if(partitions[i][j].crosses.size()>2){
//							local[k].multiple_crosses_count++;
//						}
//					}
//				}
//				poly->reset_grid_partition();
//			}
//			targets.clear();
//			distance_targets.clear();
//		}
//		delete elem;
//	}
//
//	pthread_mutex_lock(&report_lock);
//	for(int i=0;i<totalcount;i++){
//		ctx->infos[i].merge(local[i]);
//	}
//	ctx->total_num_vertices += num_vertices;
//
//	pthread_mutex_unlock(&report_lock);
//
//	return NULL;
//}
//
//
//
//int main(int argc, char** argv) {
//	string source_path;
//	int num_threads = get_num_threads();
//	partition_context ctx;
//	int max_vpr = 100;
//	int max_polygons = INT_MAX;
//	po::options_description desc("query usage");
//	desc.add_options()
//		("help,h", "produce help message")
//		("num_threads,t", po::value<int>(&num_threads), "number of threads")
//		("min_vpr", po::value<int>(&ctx.vpr), "minimum vpr")
//		("max_vpr", po::value<int>(&max_vpr), "maximum vpr")
//		("max_polygons,p", po::value<int>(&max_polygons), "maximum number of polygons")
//		("query_count,q", po::value<int>(&num_queries), "number of queries for each polygon")
//		("source_path,s", po::value<string>(&source_path), "the path to the source")
//		;
//	po::variables_map vm;
//	po::store(po::parse_command_line(argc, argv, desc), vm);
//	if (vm.count("help")) {
//		cout << desc << "\n";
//		return 0;
//	}
//	po::notify(vm);
//	for(int i=ctx.vpr;i<=max_vpr;i++){
//		ctx.vprs.push_back(i);
//		ctx.infos.push_back(global_info());
//	}
//
//	global_space = MyPolygon::gen_box(-180,-90,180,90);
//	timeval start = get_cur_time();
//
//	pthread_t threads[num_threads];
//
//	for(int i=0;i<num_threads;i++){
//		pthread_create(&threads[i], NULL, partition, (void *)&ctx);
//	}
//
//	ifstream is;
//	is.open(source_path.c_str(), ios::in | ios::binary);
//	int eid = 0;
//	queue_element *elem = new queue_element(eid++);
//	int loaded = 0;
//	size_t pid = 0;
//	while(!is.eof()&&max_polygons-->0){
//		MyPolygon *poly = MyPolygon::read_polygon_binary_file(is);
//		if(!poly||poly->get_num_vertices()<500){
//			if(poly){
//				delete poly;
//			}
//			continue;
//		}
//		poly->setid(pid++);
//		elem->polys.push_back(poly);
//		if(elem->polys.size()==element_size){
//			while(shared_queue.size()>10*num_threads){
//				usleep(20);
//			}
//			pthread_mutex_lock(&poly_lock);
//			shared_queue.push(elem);
//			pthread_mutex_unlock(&poly_lock);
//			elem = new queue_element(eid++);
//		}
//		if(++loaded%10000==0){
//			log("processed %d polygons",loaded);
//		}
//	}
//	if(elem->polys.size()>0){
//		pthread_mutex_lock(&poly_lock);
//		shared_queue.push(elem);
//		pthread_mutex_unlock(&poly_lock);
//	}else{
//		delete elem;
//	}
//	stop = true;
//
//	for(int i = 0; i < num_threads; i++ ){
//		void *status;
//		pthread_join(threads[i], &status);
//	}
//	logt("processed %d polygons",start,loaded);
//	for(int i=0;i<ctx.vprs.size();i++){
//		global_info info = ctx.infos[i];
//		int total = info.ins+info.outs+info.borders;
//		assert(total>0);
//
//		printf("%d\t%ld\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%f\n",
//				ctx.vprs[i],
//				info.num_pixels/pid,
//				1.0*info.borders/total,
//				(0.5*info.crosses_count*8+4.25*info.num_pixels)/(ctx.total_num_vertices*16),
//				0.5*info.crosses_count/info.line_count,
//				1.0*info.crosses_count/info.borders,
//				info.partition_time/pid,
//				info.query_time/pid,
//				info.nopartition_query_time/pid,
//				info.distance_query_time/pid,
//				info.nopartition_distance_query_time/pid,
//				info.two_crosses_count,
//				info.multiple_crosses_count,
//				info.two_crosses_count*1.0/(info.two_crosses_count+info.multiple_crosses_count)
//		);
//	}
//
//	delete global_space;
//
//	return 0;
//}
//
//
//
