/*
 * Parser.cpp
 *
 *  Created on: May 9, 2020
 *      Author: teng
 */



#include "../index/RTree.h"
#include <queue>
#include <fstream>
#include "../geometry/MyPolygon.h"


class polygon_list{
public:
	int id = 0;
	vector<MyPolygon *> polygons;
	pthread_mutex_t lk;
	polygon_list(){
		pthread_mutex_init(&lk, NULL);
	}

};
// some shared parameters

RTree<polygon_list *, double, 2, double> tree;

bool MySearchCallback(polygon_list *polys, void* arg){
	query_context *ctx = (query_context *)arg;

	MyPolygon *target = (MyPolygon *)ctx->target;
	pthread_mutex_lock(&polys->lk);
    polys->polygons.push_back(target);
	pthread_mutex_unlock(&polys->lk);

	return true;
}

void *query(void *args){
	query_context *ctx = (query_context *)args;
	query_context *gctx = ctx->global_ctx;
	//log("thread %d is started",ctx->thread_id);
	while(ctx->next_batch(100)){
		for(int i=ctx->index;i<ctx->index_end;i++){
			//gctx->source_polygons[i]->getMBB()->print();
			ctx->target = (void *)gctx->source_polygons[i];
			tree.Search(gctx->source_polygons[i]->getMBB()->low, gctx->source_polygons[i]->getMBB()->high, MySearchCallback, (void *)ctx);
			ctx->report_progress();
		}
	}
	ctx->merge_global();

	return NULL;
}


vector<MyPolygon *> local_load_binary_file(const char *path, query_context &ctx){
	vector<MyPolygon *> polygons;
	if(!file_exist(path)){
		log("%s does not exist",path);
		return polygons;
	}
	log("loading polygon from %s",path);
	struct timeval start = get_cur_time();
	ifstream infile;
	infile.open(path, ios::in | ios::binary);
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
	num_polygons = min(num_polygons, ctx.max_num_polygons);

	log("start loading %d polygons",num_polygons);

	size_t num_edges = 0;
	size_t data_size = 0;
	size_t next = 10;
	size_t id = 0;
	for(size_t i=0;i<num_polygons;i++){
		if(i*100/num_polygons>=next){
			log("loaded %d%% %ld",next,id);
			next+=10;
		}

		infile.seekg(offsets[i], infile.beg);
		MyPolygon *poly = MyPolygon::read_polygon_binary_file(infile);
		assert(poly);

		if(poly->get_num_vertices()<100 && !tryluck(ctx.sample_rate)){
			delete poly;
			continue;
		}

		num_edges += poly->get_num_vertices();
		data_size += poly->get_data_size();
		poly->setid(id++);
		polygons.push_back(poly);

	}
	infile.close();

	logt("loaded %ld polygons each with %ld edges %ld MB", start, polygons.size(),num_edges/polygons.size(),data_size/1024/1024);
	return polygons;
}

int main(int argc, char** argv) {

	MyMultiPolygon *mp = MyMultiPolygon::read_multipolygon();
	vector<MyPolygon *> polygons = mp->get_polygons();

	vector<polygon_list *> partitions;
	for(int i=0;i<polygons.size();i++){
		polygon_list *pos = new polygon_list();
		pos->id = i;
		partitions.push_back(pos);
		polygons[i]->setid(i);
		tree.Insert(polygons[i]->getMBB()->low, polygons[i]->getMBB()->high, partitions[i]);
	}

	query_context global_ctx;
	global_ctx = get_parameters(argc, argv);
	global_ctx.source_polygons = local_load_binary_file(global_ctx.source_path.c_str(),global_ctx);
	global_ctx.target_num = global_ctx.source_polygons.size();
	timeval start = get_cur_time();

	// read all the points
	start = get_cur_time();
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

	for(int i = 0; i < global_ctx.num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}
	logt("total query",start);

	int total = 0;

	int max_one = 0;
	for(int i=0;i<partitions.size();i++){
		char path[256];
		sprintf(path, "part_polygons/%d.dat",i);

		dump_polygons_to_file(partitions[i]->polygons, path);

		log("%d %d",i,partitions[i]->polygons.size());
		if(partitions[i]->polygons.size()>partitions[max_one]->polygons.size()){
			max_one = i;
		}
	}
	log("%d %d",max_one,partitions[max_one]->polygons.size());

//	for(MyPolygon *p:partitions[max_one]->polygons){
//		p->print(false, false);
//	}
	//polygons[max_one]->print(false, false);

	return 0;
}



