#include "../geometry/MyPolygon.h"
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <vector>
#include <string.h>
#include <fstream>
#include <iostream>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

using namespace std;

#define MAX_THREAD_NUM 100

// some shared parameters
string processing_line[MAX_THREAD_NUM];
bool is_working[MAX_THREAD_NUM];
pthread_mutex_t output_lock;
bool stop = false;
int big_threshold = 100;
int small_threshold = 1000;

ofstream big_wf;
ofstream small_wf;

void *process_wkt(void *args){

	int id = *(int *)args;
	pthread_mutex_lock(&output_lock);
	log("thread %d is started", id);
	pthread_mutex_unlock(&output_lock);
	vector<MyMultiPolygon *> mpolys;
	size_t data_size_big = 0;
	size_t data_size_small = 0;
	size_t buffer_size = 100*1024*1024;
	char *data_buffer_big = new char[buffer_size];
	char *data_buffer_small = new char[buffer_size];

	//Parsing polygon input
	while (!stop||is_working[id]) {
		if(!is_working[id]){
			usleep(10);
			continue;
		}
		MyMultiPolygon *mp = new MyMultiPolygon(processing_line[id].c_str());

		vector<MyPolygon *> polygons = mp->get_polygons();
		for(MyPolygon *p:polygons){
			if(p->get_num_vertices()>=big_threshold){
				if(p->get_data_size()+data_size_big>buffer_size){
					pthread_mutex_lock(&output_lock);
					big_wf.write(data_buffer_big, data_size_big);
					pthread_mutex_unlock(&output_lock);
					data_size_big = 0;
				}
				data_size_big += p->encode_to(data_buffer_big+data_size_big);
			}
			if(p->get_num_vertices()<=small_threshold){
				if(p->get_data_size()+data_size_small>buffer_size){
					pthread_mutex_lock(&output_lock);
					small_wf.write(data_buffer_small, data_size_small);
					pthread_mutex_unlock(&output_lock);
					data_size_small = 0;
				}
				data_size_small += p->encode_to(data_buffer_small+data_size_small);
			}
		}
		delete mp;
		processing_line[id].clear();
		is_working[id] = false;
	} // end of while
	pthread_mutex_lock(&output_lock);
	big_wf.write(data_buffer_big, data_size_big);
	small_wf.write(data_buffer_small, data_size_small);
	pthread_mutex_unlock(&output_lock);
	delete []data_buffer_big;
	delete []data_buffer_small;

	pthread_exit(NULL);
	return NULL;
}


int main(int argc, char** argv) {
	string path1;
	string path2;

	po::options_description desc("load usage");
	desc.add_options()
		("help,h", "produce help message")
		("path_big", po::value<string>(&path1)->required(), "path for the big polygons")
		("path_small", po::value<string>(&path2)->required(), "path for the small polygons")
		("big_threshold,b", po::value<int>(&big_threshold), "minimum number of vertices for big polygons")
		("small_threshold,s", po::value<int>(&small_threshold), "maximum number of vertices for small polygons")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
		cout << desc << "\n";
		return 0;
	}
	po::notify(vm);

	big_wf.open(path1.c_str(), ios::out | ios::binary |ios::trunc);
	small_wf.open(path2.c_str(), ios::out | ios::binary|ios::trunc);

	int num_threads = get_num_threads();
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

	while(getline(cin, input_line)) {
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
	big_wf.close();
	small_wf.close();
}
