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

vector<size_t> big_offsets;
vector<size_t> small_offsets;

size_t cur_big_offset = 0;
size_t cur_small_offset = 0;
bool fix = false;

void *process_wkt(void *args){

	int id = *(int *)args;
	pthread_mutex_lock(&output_lock);
	log("thread %d is started", id);
	pthread_mutex_unlock(&output_lock);
	size_t data_size_big = 0;
	size_t data_size_small = 0;
	size_t buffer_size = 100*1024*1024;
	char *data_buffer_big = new char[buffer_size];
	char *data_buffer_small = new char[buffer_size];

	vector<size_t> big_size;
	vector<size_t> small_size;;

	//Parsing polygon input
	while (!stop||is_working[id]) {
		if(!is_working[id]){
			usleep(10);
			continue;
		}
		MyMultiPolygon *mp = new MyMultiPolygon(processing_line[id].c_str());

		vector<MyPolygon *> polygons = mp->get_polygons();
		for(MyPolygon *p:polygons){
			if(p->get_num_vertices()>=big_threshold || p->get_num_vertices()<=small_threshold){
				if(fix){
					p->boundary->fix();
				}
			}
			//log("processed polygon with %d vertices", p->get_num_vertices());
			if(p->get_num_vertices()>=big_threshold){
				if(p->get_data_size()+data_size_big>buffer_size){
					pthread_mutex_lock(&output_lock);
					big_wf.write(data_buffer_big, data_size_big);
					for(size_t s:big_size){
						big_offsets.push_back(cur_big_offset);
						cur_big_offset += s;
					}
					pthread_mutex_unlock(&output_lock);
					data_size_big = 0;
					big_size.clear();
				}
				big_size.push_back(p->get_data_size());
				data_size_big += p->encode_to(data_buffer_big+data_size_big);
			}
			if(p->get_num_vertices()<=small_threshold){
				if(p->get_data_size()+data_size_small>buffer_size){
					pthread_mutex_lock(&output_lock);
					small_wf.write(data_buffer_small, data_size_small);
					for(size_t s:small_size){
						small_offsets.push_back(cur_small_offset);
						cur_small_offset += s;
					}
					pthread_mutex_unlock(&output_lock);
					data_size_small = 0;
					small_size.clear();
				}
				small_size.push_back(p->get_data_size());
				data_size_small += p->encode_to(data_buffer_small+data_size_small);
			}
		}
		delete mp;
		processing_line[id].clear();
		is_working[id] = false;
	} // end of while
	pthread_mutex_lock(&output_lock);
	if(data_size_big>0){
		big_wf.write(data_buffer_big, data_size_big);
		for(size_t s:big_size){
			big_offsets.push_back(cur_big_offset);
			cur_big_offset += s;
		}
		big_size.clear();
	}
	if(data_size_small>0){
		small_wf.write(data_buffer_small, data_size_small);
		for(size_t s:small_size){
			small_offsets.push_back(cur_small_offset);
			cur_small_offset += s;
		}
		small_size.clear();
	}

	pthread_mutex_unlock(&output_lock);
	delete []data_buffer_big;
	delete []data_buffer_small;

	pthread_exit(NULL);
	return NULL;
}


int main(int argc, char** argv) {
	string path1;
	string path2;
	int num_threads = get_num_threads();

	po::options_description desc("load usage");
	desc.add_options()
		("help,h", "produce help message")
		("path_big", po::value<string>(&path1)->required(), "path for the big polygons")
		("path_small", po::value<string>(&path2)->required(), "path for the small polygons")
		("big_threshold,b", po::value<int>(&big_threshold), "minimum number of vertices for big polygons")
		("small_threshold,s", po::value<int>(&small_threshold), "maximum number of vertices for small polygons")
		("num_threads,t", po::value<int>(&num_threads), "number of threads")
		("fix,f", "fix the boundary while loading")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
		cout << desc << "\n";
		return 0;
	}
	po::notify(vm);
	fix = vm.count("fix");

	big_wf.open(path1.c_str(), ios::out | ios::binary |ios::trunc);
	small_wf.open(path2.c_str(), ios::out | ios::binary|ios::trunc);

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
		if(num_objects%100000==0){
			log("processed %d objects", num_objects);
		}
	}
	stop = true;

	for(int i = 0; i < num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}

	big_wf.write((char *)&big_offsets[0], sizeof(size_t)*big_offsets.size());
	size_t bs = big_offsets.size();
	big_wf.write((char *)&bs, sizeof(size_t));

	small_wf.write((char *)&small_offsets[0], sizeof(size_t)*small_offsets.size());
	size_t ss = small_offsets.size();
	small_wf.write((char *)&ss, sizeof(size_t));

	logt("processed %d objects %ld big and %ld small", start_time, num_objects, bs, ss);
	big_wf.close();
	small_wf.close();
}
