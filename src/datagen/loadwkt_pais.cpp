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
queue<pair<string, string>> input_lines;
bool stop = false;

pthread_mutex_t output_lock;
pthread_mutex_t fetch_lock;

class ImageMeta{
public:
	string path;
	ofstream os;
	vector<PolygonMeta> pmeta;
	size_t global_offset = 0;
	string image_id;
	pthread_mutex_t lk;

	void lock(){
		pthread_mutex_lock(&lk);
	}
	void unlock(){
		pthread_mutex_unlock(&lk);
	}
	ImageMeta(const char *pt){
		path = pt;
		os.open(pt, ios::out | ios::binary |ios::trunc);
		assert(os.is_open());
		pthread_mutex_init(&lk, NULL);
	}
	~ImageMeta(){
		pmeta.clear();
		os.close();
	}
	void persist(MyPolygon *poly, char *buffer){
		size_t s = poly->encode(buffer);
		PolygonMeta pm = poly->get_meta();

		lock();
		pm.offset = global_offset;
		global_offset += s;
		os.write(buffer, s);
		pmeta.push_back(pm);
		unlock();
	}
	void finish(){
		os.write((char *)&pmeta[0], sizeof(PolygonMeta)*pmeta.size());
		size_t bs = pmeta.size();
		os.write((char *)&bs, sizeof(size_t));
	}

};

map<string, ImageMeta *> images;


bool fix = false;

size_t valid_line_count = 0;

void *process_wkt(void *args){

	int id = *(int *)args;
	pthread_mutex_lock(&output_lock);
	log("thread %d is started", id);
	pthread_mutex_unlock(&output_lock);

	char *buffer = new char[buffer_size];
	size_t valid_local = 0;
	//Parsing polygon input
	while (!stop) {
		pair<string, string> linepair;
		bool isempty = true;
		pthread_mutex_lock(&fetch_lock);
		if(!input_lines.empty()){
			linepair = input_lines.front();
			input_lines.pop();
			isempty = false;
		}
		pthread_mutex_unlock(&fetch_lock);
		if(isempty){
			usleep(1);
			continue;
		}
		string wkt = linepair.second;
		string image_id = linepair.first;

		if(MyMultiPolygon::validate_wkt(wkt)){
			MyMultiPolygon *mp = new MyMultiPolygon(wkt.c_str());
			vector<MyPolygon *> polygons = mp->get_polygons();
			for(MyPolygon *p:polygons){
				if(fix){
					p->boundary->fix();
				}
				images[image_id]->persist(p, buffer);
				//log("%f",p->getMBB()->area());
			}
			delete mp;
			valid_local++;
		}else{
			printf("%s\n",wkt.c_str());
			log("invalid wkt");
		}
		wkt.clear();
		image_id.clear();
	} // end of while
	pthread_mutex_lock(&output_lock);
	valid_line_count += valid_local;
	pthread_mutex_unlock(&output_lock);

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

	char path[256];
	while(getline(is, input_line)) {
		vector<string> fields;
		tokenize(input_line, fields, ",", true);

		if(fields.size()<3){
			fields.clear();
			continue;
		}
		if(images.find(fields[0])==images.end()){
			sprintf(path,"%s/%s.idl",outpath.c_str(),fields[0].c_str());
			ImageMeta *img = new ImageMeta(path);
			images[fields[0]] = img;
		}

		pthread_mutex_lock(&fetch_lock);
		input_lines.push(pair<string, string>(fields[0], fields[2]));
		pthread_mutex_unlock(&fetch_lock);
		if(++num_objects%100==0){
			log_refresh("processed %d objects", num_objects);
		}
	}
	stop = true;

	for(int i = 0; i < num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}

	for(auto it = images.begin();it!=images.end();it++){
		it->second->finish();
		delete it->second;
	}
	images.clear();

	logt("processed %d lines %ld valid %ld invalid", start_time, num_objects, valid_line_count, num_objects - valid_line_count);
}
