#include "MyPolygon.h"
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <vector>
#include <string.h>
/* 
 * RESQUE processing engine v3.0
 *   It supports spatial join and nearest neighbor query with different predicates
 *   1) parseParameters
 *   2) readCacheFile - metadata such as partition schemata
 *   3) for every input line in the current tile
 *         an input line represents an object
 *         save geometry and original data in memory
 *         execute join operation when finish reading a tile
 *   4) Join operation between 2 sets or a single set
 *         build Rtree index on the second data set
 *         for every object in the first data set
 *            using Rtree index of the second data set
 *              check for MBR/envelope intersection
 *              output the pair result or save pair statistics
 *   5) Output final statistics (opt)
 *   Requirement (input files): see the Wiki
 * */

using namespace std;

#define QUEUE_SIZE 100
#define MAX_THREAD_NUM 100

// some shared parameters
string processing_line[MAX_THREAD_NUM];
bool is_working[MAX_THREAD_NUM];
pthread_mutex_t line_lock;
pthread_mutex_t output_lock;
bool stop = false;

MyPolygon *max_poly = NULL;
long total_num_vertices = 0;
long total_num_polygons = 0;

long *vertices_count;
int num_stats = 50000;


void *process_wkt(void *args){

	long *local_vertices_count = new long[num_stats];
	int id = *(int *)args;
	pthread_mutex_lock(&output_lock);
	log("thread %d is started", id);
	pthread_mutex_unlock(&output_lock);
	MyPolygon *local_max_poly = NULL;
	long local_num_vertices = 0;
	long local_num_polygons = 0;
	long local_id = 0;

	// Reading from the cache file
	/* Parsing polygon input */
	while (!stop||is_working[id]) {
		if(!is_working[id]){
			usleep(10);
			continue;
		}
		MyMultiPolygon *mp = new MyMultiPolygon(processing_line[id].c_str());
		vector<MyPolygon *> polygons = mp->get_polygons();
		for(MyPolygon *p:polygons){
			local_num_polygons++;
			local_num_vertices += p->num_boundary_vertices();
			if(!local_max_poly){
				local_max_poly = p->clone();
			}else if(p->num_boundary_vertices()>local_max_poly->num_boundary_vertices()){
				delete local_max_poly;
				local_max_poly = p->clone();
			}
			if(p->num_boundary_vertices()>num_stats){
				local_vertices_count[num_stats-1]++;
			}else{
				local_vertices_count[p->num_boundary_vertices()-1]++;
			}
		}
		delete mp;
		processing_line[id].clear();
		is_working[id] = false;
	} // end of while
	if(local_max_poly){
		pthread_mutex_lock(&output_lock);
		if(!max_poly){
			max_poly = local_max_poly->clone();
		}else if(local_max_poly->num_boundary_vertices()>max_poly->num_boundary_vertices()){
			delete max_poly;
			max_poly = local_max_poly;
		}
		total_num_polygons += local_num_polygons;
		total_num_vertices += local_num_vertices;
		for(int i=0;i<num_stats;i++){
			vertices_count[i] += local_vertices_count[i];
		}
		delete []local_vertices_count;
		pthread_mutex_unlock(&output_lock);
	}
	pthread_exit(NULL);
	return NULL;
}


int main(int argc, char** argv) {
	int num_threads = get_num_threads();
	if(argc>=2){
		num_threads = atoi(argv[1]);
		if(num_threads == 0){
			num_threads = get_num_threads();
		}
	}
	if(argc>=3){
		num_stats = atoi(argv[2]);
	}
	vertices_count = new long[num_stats];


	pthread_t threads[num_threads];
	int id[num_threads];
	for(int i=0;i<num_threads;i++){
		id[i] = i;
		processing_line[i].clear();
		is_working[i] = false;
	}
	for(int i=0;i<num_threads;i++){
		pthread_create(&threads[i], NULL, process_wkt, (void *)&id[i]);
	}
	struct timeval start_time = get_cur_time();
	std::string input_line;
	int num_objects = 0;
	int next_report = 0;
	long processed_size = 0;

	while (getline(cin, input_line)) {
		while(true){
			bool assigned = false;
			for(int i=0;i<num_threads;i++){
				if(is_working[i]==false){
					processing_line[i] = input_line;
					is_working[i] = true;
					assigned = true;
					break;
				}
			}
			if(assigned){
				break;
			}
			usleep(10);
		}
		processed_size += input_line.size()+1;
		num_objects++;
		if(num_objects%1000000==0){
			log("processed %d objects", num_objects);
		}
	}
	stop = true;

	for(int i = 0; i < num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}

	logt("processed %d objects", start_time, num_objects);


	if(max_poly){
		//max_poly->print();
		cout<<max_poly->num_boundary_vertices()<<endl;
		delete max_poly;
	}
	cerr<<total_num_polygons<<endl;
	cerr<<total_num_vertices<<endl;
	cerr<<total_num_vertices/total_num_polygons<<endl;

	long cum = 0;
	for(int i=0;i<num_stats;i++){
		if(vertices_count[i]!=0){
			cum += vertices_count[i];
			cout<<i<<","<<vertices_count[i]<<","<<(double)cum/total_num_polygons<<endl;
		}
	}
	delete []vertices_count;

	pthread_exit(NULL);
	return true;
}
