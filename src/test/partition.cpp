/*
 * Parser.cpp
 *
 *  Created on: May 9, 2020
 *      Author: teng
 */



#include "../index/RTree.h"
#include <queue>
#include <fstream>
#include <iostream>
#include "../include/MyPolygon.h"
#include <queue>
using namespace std;


queue<MyPolygon *> poly_queue;
bool finished = false;
pthread_mutex_t qlock;
const char *target_folder;

class polygon_list{
public:
	int id = 0;
	vector<size_t> offsets;
	size_t cur_offset = 0;
	ofstream fs;
	pthread_mutex_t lk;

	void push(MyPolygon *poly, char *buffer){
		size_t size = poly->encode_to(buffer);
		pthread_mutex_lock(&lk);
		offsets.push_back(cur_offset);
		cur_offset += size;
		fs.write(buffer, size);
		pthread_mutex_unlock(&lk);
	}

	polygon_list(int i){
		id = i;
		char path[256];
		sprintf(path, "%s/%d.dat",target_folder, i);
		fs.open(path, ios::out | ios::binary |ios::trunc);
		assert(fs.is_open());
		pthread_mutex_init(&lk, NULL);
	}

	~polygon_list(){
		size_t ss = offsets.size();
		if(ss>0){
			fs.write((char *)&offsets[0], sizeof(size_t)*ss);
		}
		fs.write((char *)&ss, sizeof(size_t));
		fs.close();
	}

};
// some shared parameters
RTree<polygon_list *, double, 2, double> tree;


bool MySearchCallback(polygon_list *polys, void* arg){
	query_context *ctx = (query_context *)arg;

	MyPolygon *target = (MyPolygon *)ctx->target;
	char *buffer = (char *)ctx->target2;
	polys->push(target, buffer);

	return true;
}


void *query(void *args){
	char *buffer = new char[100*1024*1024];
	query_context *ctx = (query_context *)args;
	ctx->target2 = (void *)buffer;
	//log("thread %d is started",ctx->thread_id);
	while(!finished){
		MyPolygon *poly = NULL;
		pthread_mutex_lock(&qlock);
		if(poly_queue.size()>0){
			poly = poly_queue.front();
			poly_queue.pop();
		}
		pthread_mutex_unlock(&qlock);
		if(poly!=NULL){
			ctx->target = (void *)poly;
			tree.Search(poly->getMBB()->low, poly->getMBB()->high, MySearchCallback, (void *)ctx);
			delete poly;
		}else{
			usleep(10);
		}
	}
	return NULL;
}

int main(int argc, char** argv) {

	query_context global_ctx;
	global_ctx = get_parameters(argc, argv);
	target_folder = global_ctx.target_path.c_str();
	// generate the index with the input partitions
	MyMultiPolygon *mp = MyMultiPolygon::read_multipolygon();
	vector<MyPolygon *> polygons = mp->get_polygons();
	vector<polygon_list *> partitions;
	for(int i=0;i<polygons.size();i++){
		polygon_list *pos = new polygon_list(i);
		partitions.push_back(pos);
		polygons[i]->setid(i);
		tree.Insert(polygons[i]->getMBB()->low, polygons[i]->getMBB()->high, partitions[i]);
	}

	timeval start = get_cur_time();

	pthread_t threads[global_ctx.num_threads];
	query_context ctx[global_ctx.num_threads];
	for(int i=0;i<global_ctx.num_threads;i++){
		ctx[i] = global_ctx;
		ctx[i].thread_id = i;
		ctx[i].global_ctx = &global_ctx;
	}
	for(int i=0;i<global_ctx.num_threads;i++){
		pthread_create(&threads[i], NULL, query, (void *)&ctx[i]);
	}

	// the main thread load the target
	{
		if(!file_exist(global_ctx.source_path.c_str())){
			log("%s does not exist",global_ctx.source_path.c_str());
			return 0;
		}
		log("loading polygon from %s",global_ctx.source_path.c_str());
		ifstream infile;
		infile.open(global_ctx.source_path.c_str(), ios::in | ios::binary);
		size_t off;
		size_t num_polygons = 0;
		infile.seekg(0, infile.end);
		//seek to the first polygon
		infile.seekg(-sizeof(size_t), infile.end);
		infile.read((char *)&num_polygons, sizeof(size_t));
		assert(num_polygons>0 && "the file should contain at least one polygon");
		size_t *offsets = new size_t[num_polygons];

		infile.seekg(-sizeof(size_t)*(num_polygons+1), infile.end);
		infile.read((char *)offsets, sizeof(size_t)*num_polygons);
		num_polygons = min(num_polygons, global_ctx.max_num_polygons);

		log("start loading %d polygons",num_polygons);

		size_t num_edges = 0;
		size_t data_size = 0;
		size_t next = 1;
		for(size_t i=0;i<num_polygons;i++){
			if(i*100/num_polygons>=next){
				logt("loaded %d(%d%%)",start, i, next);
				next+=1;
			}

			infile.seekg(offsets[i], infile.beg);
			MyPolygon *poly = MyPolygon::read_polygon_binary_file(infile);
			poly->setid(i);
			assert(poly);
			pthread_mutex_lock(&qlock);
			poly_queue.push(poly);
			pthread_mutex_unlock(&qlock);

			num_edges += poly->get_num_vertices();
			data_size += poly->get_data_size();
		}
		infile.close();
		finished = true;

		logt("loaded %ld polygons each with %ld edges %ld MB", start, num_polygons,num_edges/num_polygons,data_size/1024/1024);
	}

	for(int i = 0; i < global_ctx.num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}
	logt("total query",start);

	int total = 0;

	int max_one = 0;
	int max_size = 0;
	for(int i=0;i<partitions.size();i++){
		char path[256];
		sprintf(path, "part_polygons/%d.dat",i);
		size_t s = partitions[i]->offsets.size();
		delete partitions[i];
		log("%d %d",i,s);
		if(max_size<s){
			max_one = i;
			max_size = s;
		}
	}
	log("%d %d",max_one,max_size);

//	for(MyPolygon *p:partitions[max_one]->polygons){
//		p->print(false, false);
//	}
	//polygons[max_one]->print(false, false);

	return 0;
}



