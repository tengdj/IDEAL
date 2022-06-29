#include "../include/MyPolygon.h"
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <vector>
#include <string.h>
#include <fstream>
#include <iostream>
#include <boost/program_options.hpp>
#include <queue>

namespace po = boost::program_options;

using namespace std;

#define MAX_THREAD_NUM 100

const size_t buffer_size = 100*1024*1024;

// some shared parameters
queue<string> input_lines;
bool stop = false;

pthread_mutex_t output_lock;
pthread_mutex_t fetch_lock;

ofstream os;
vector<PolygonMeta> pmeta;

size_t global_offset = 0;
bool fix = false;

size_t valid_line_count = 0;

void *process_wkt(void *args){

	int id = *(int *)args;
	pthread_mutex_lock(&output_lock);
	log("thread %d is started", id);
	pthread_mutex_unlock(&output_lock);
	size_t data_size = 0;
	char *data_buffer = new char[buffer_size];

	vector<PolygonMeta> local_pmeta;
	size_t valid_local = 0;
	//Parsing polygon input
	while (!stop) {
		string line;
		bool isempty = true;
		pthread_mutex_lock(&fetch_lock);
		if(!input_lines.empty()){
			line = input_lines.front();
			input_lines.pop();
			isempty = false;
		}
		pthread_mutex_unlock(&fetch_lock);
		if(isempty){
			usleep(1);
			continue;
		}
		if(MyMultiPolygon::validate_wkt(line)){
			MyMultiPolygon *mp = new MyMultiPolygon(line.c_str());

			vector<MyPolygon *> polygons = mp->get_polygons();
			for(MyPolygon *p:polygons){
				if(fix){
					p->boundary->fix();
				}
				//log("processed polygon with %d vertices", p->get_num_vertices());
				if(p->get_data_size() + data_size>buffer_size){
					pthread_mutex_lock(&output_lock);
					os.write(data_buffer, data_size);
					for(PolygonMeta &m:local_pmeta){
						m.offset = global_offset;
						global_offset += m.size;
					}
					pmeta.insert(pmeta.end(), local_pmeta.begin(), local_pmeta.end());
					pthread_mutex_unlock(&output_lock);
					data_size = 0;
					local_pmeta.clear();
				}
				data_size += p->encode(data_buffer + data_size);
				local_pmeta.push_back(p->get_meta());
			}
			delete mp;
			valid_local++;
		}else{
			printf("%s\n",line.c_str());
			log("invalid wkt");
		}
		line.clear();
	} // end of while
	if(data_size>0){
		pthread_mutex_lock(&output_lock);
		os.write(data_buffer, data_size);
		for(PolygonMeta &m:local_pmeta){
			m.offset = global_offset;
			global_offset += m.size;
		}
		pmeta.insert(pmeta.end(), local_pmeta.begin(), local_pmeta.end());
		pthread_mutex_unlock(&output_lock);
	}
	pthread_mutex_lock(&output_lock);
	valid_line_count += valid_local;
	pthread_mutex_unlock(&output_lock);

	delete []data_buffer;

	pthread_exit(NULL);
	return NULL;
}


int main(int argc, char** argv) {
	string inpath;
	string outpath;
	int num_threads = get_num_threads();

	po::options_description desc("load usage");
	desc.add_options()
		("help,h", "produce help message")
		("input,i", po::value<string>(&inpath)->required(), "path for the big polygons")
		("output,o", po::value<string>(&outpath)->required(), "path for the small polygons")
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

	os.open(outpath.c_str(), ios::out | ios::binary |ios::trunc);
	assert(os.is_open());
	std::ifstream is(inpath.c_str());
	assert(is.is_open());

	pthread_t threads[num_threads];
	int id[num_threads];
	for(int i=0;i<num_threads;i++){
		id[i] = i;
	}
	for(int i=0;i<num_threads;i++){
		pthread_create(&threads[i], NULL, process_wkt, (void *)&id[i]);
	}
	struct timeval start_time = get_cur_time();
	std::string input_line;
	size_t num_objects = 0;
	size_t processed_size = 0;

	while(getline(is, input_line)) {
		pthread_mutex_lock(&fetch_lock);
		input_lines.push(input_line);
		pthread_mutex_unlock(&fetch_lock);
		processed_size += input_line.size()+1;
		if(++num_objects%100==0){
			log_refresh("processed %d objects", num_objects);
		}
	}
	stop = true;

	for(int i = 0; i < num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}

	os.write((char *)&pmeta[0], sizeof(PolygonMeta)*pmeta.size());
	size_t bs = pmeta.size();
	os.write((char *)&bs, sizeof(size_t));

	logt("processed %d lines %ld valid %ld invalid", start_time, num_objects, valid_line_count, num_objects - valid_line_count);
	os.close();
	pmeta.clear();
}
