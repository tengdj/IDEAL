//#include "../include/MyPolygon.h"
//#include <stdio.h>
//#include <stdlib.h>
//#include <pthread.h>
//#include <vector>
//#include <string.h>
///*
// *
// * collect the statistic of the input data
// *
// * */
//
//using namespace std;
//
//#define QUEUE_SIZE 100
//#define MAX_THREAD_NUM 100
//
//// some shared parameters
//string processing_line[MAX_THREAD_NUM];
//bool is_working[MAX_THREAD_NUM];
//pthread_mutex_t line_lock;
//pthread_mutex_t output_lock;
//bool stop = false;
//
//vector<MyPolygon *> large_polys;
//int large_poly_threshold = 1000;
//MyPolygon *max_poly;
//int max_num_polygs = 100;
//long total_num_vertices = 0;
//long total_num_polygons = 0;
//
//long *vertices_count;
//int num_stats = 50000;
//
//
//void *process_wkt(void *args){
//
//	long *local_vertices_count = new long[num_stats];
//	int id = *(int *)args;
//	pthread_mutex_lock(&output_lock);
//	log("thread %d is started", id);
//	pthread_mutex_unlock(&output_lock);
//	MyPolygon *local_max_poly = NULL;
//	long local_num_vertices = 0;
//	long local_num_polygons = 0;
//	long local_id = 0;
//
//	// Reading from the cache file
//	/* Parsing polygon input */
//	while (!stop||is_working[id]) {
//		if(!is_working[id]){
//			usleep(10);
//			continue;
//		}
//		MyMultiPolygon *mp = new MyMultiPolygon(processing_line[id].c_str());
//		vector<MyPolygon *> polygons = mp->get_polygons();
//		for(MyPolygon *p:polygons){
//			local_num_polygons++;
//			local_num_vertices += p->get_num_vertices();
//			if(!local_max_poly){
//				local_max_poly = p->clone();
//			}else if(p->get_num_vertices()>local_max_poly->get_num_vertices()){
//				delete local_max_poly;
//				local_max_poly = p->clone();
//			}
//			if(p->get_num_vertices()>num_stats){
//				local_vertices_count[num_stats-1]++;
//			}else{
//				local_vertices_count[p->get_num_vertices()-1]++;
//			}
//			if(p->get_num_vertices()>large_poly_threshold){
//				pthread_mutex_lock(&output_lock);
//				MyPolygon *tp = p->clone();
//				large_polys.push_back(tp);
//				pthread_mutex_unlock(&output_lock);
//			}
//		}
//		delete mp;
//		processing_line[id].clear();
//		is_working[id] = false;
//	} // end of while
//	if(local_max_poly){
//		pthread_mutex_lock(&output_lock);
//		if(!max_poly){
//			max_poly = local_max_poly->clone();
//		}else if(local_max_poly->get_num_vertices()>max_poly->get_num_vertices()){
//			delete max_poly;
//			max_poly = local_max_poly;
//		}
//		total_num_polygons += local_num_polygons;
//		total_num_vertices += local_num_vertices;
//		for(int i=0;i<num_stats;i++){
//			vertices_count[i] += local_vertices_count[i];
//		}
//		delete []local_vertices_count;
//		pthread_mutex_unlock(&output_lock);
//	}
//	pthread_exit(NULL);
//	return NULL;
//}
//
//
//int collect_stats(int argc, char** argv) {
//	int num_threads = get_num_threads();
//	if(argc>=2){
//		num_threads = atoi(argv[1]);
//		if(num_threads == 0){
//			num_threads = get_num_threads();
//		}
//	}
//	if(argc>=3){
//		num_stats = atoi(argv[2]);
//	}
//	vertices_count = new long[num_stats];
//
//
//	pthread_t threads[num_threads];
//	int id[num_threads];
//	for(int i=0;i<num_threads;i++){
//		id[i] = i;
//		processing_line[i].clear();
//		is_working[i] = false;
//	}
//	for(int i=0;i<num_threads;i++){
//		pthread_create(&threads[i], NULL, process_wkt, (void *)&id[i]);
//	}
//	struct timeval start_time = get_cur_time();
//	std::string input_line;
//	int num_objects = 0;
//	int next_report = 0;
//	long processed_size = 0;
//
//	while (getline(cin, input_line)) {
//		while(true){
//			bool assigned = false;
//			for(int i=0;i<num_threads;i++){
//				if(is_working[i]==false){
//					processing_line[i] = input_line;
//					is_working[i] = true;
//					assigned = true;
//					break;
//				}
//			}
//			if(assigned){
//				break;
//			}
//			usleep(10);
//		}
//		processed_size += input_line.size()+1;
//		num_objects++;
//		if(num_objects%1000000==0){
//			log("processed %d objects", num_objects);
//		}
//	}
//	stop = true;
//
//	for(int i = 0; i < num_threads; i++ ){
//		void *status;
//		pthread_join(threads[i], &status);
//	}
//
//	logt("processed %d objects", start_time, num_objects);
//
//
//	if(max_poly){
//		//max_poly->print();
//		//cout<<max_poly->get_num_vertices()<<endl;
//		delete max_poly;
//	}
////	cerr<<total_num_polygons<<endl;
////	cerr<<total_num_vertices<<endl;
////	cerr<<total_num_vertices/total_num_polygons<<endl;
////
////	long cum = 0;
////	for(int i=0;i<num_stats;i++){
////		if(vertices_count[i]!=0){
////			cum += vertices_count[i];
////			cout<<i<<","<<vertices_count[i]<<","<<(double)cum/total_num_polygons<<endl;
////		}
////	}
//	delete []vertices_count;
//
//	return true;
//}
//
//
//int main(int argc, char **argv){
//	large_poly_threshold = 1000;
//	collect_stats(argc, argv);
//	for(int i=0;i<large_polys.size()-1;i++){
//		for(int j=large_polys.size()-1;j>i;j--){
//			if(large_polys[j]->get_num_vertices()>large_polys[j-1]->get_num_vertices()){
//				MyPolygon *tmp = large_polys[j-1];
//				large_polys[j-1] = large_polys[j];
//				large_polys[j] = tmp;
//			}
//		}
//	}
//	for(MyPolygon *p:large_polys){
//		p->print();
//		//cout<<p->get_num_vertices()<<endl;
//	}
//}
