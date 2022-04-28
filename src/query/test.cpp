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

// some shared parameters
RTree<MyPolygon *, double, 2, double> tree;

unordered_set<int> has_children;
unordered_set<int> has_parent;
bool contained = false;


bool MySearchCallback(MyPolygon *poly, void* arg){
	query_context *ctx = (query_context *)arg;
	MyPolygon *target = (MyPolygon *)ctx->target;

	if(poly->getid()==target->getid()){
		return true;
	}

	struct timeval start = get_cur_time();
	contained = false;
	if(poly->contain(target, ctx)){
		has_children.insert(poly->getid());
		has_parent.insert(target->getid());
		contained = true;
	}
	//logt("checked one %d",start,contained);

	//log("checking %d (%ld) and %d (%ld) result %d",poly->getid(),poly->get_num_vertices(),target->getid(),target->get_num_vertices(),contained);

	//ctx->found += poly->contain(*(Point *)ctx->target, ctx);

	return true;
}

void *query(void *args){
	query_context *ctx = (query_context *)args;
	query_context *gctx = ctx->global_ctx;
	//log("thread %d is started",ctx->thread_id);
	while(ctx->next_batch(100)){
		for(int i=ctx->index;i<ctx->index_end;i++){

			ctx->report_progress();
		}
	}
	ctx->merge_global();

	return NULL;
}

int main(int argc, char** argv) {

	query_context global_ctx;
	//global_ctx.max_num_polygons = 1000;
	vector<MyPolygon *> source_polygons = MyPolygon::load_binary_file("/home/teng/git/IDEAL/data/source.dat",global_ctx, false);

	timeval start = get_cur_time();
	for(MyPolygon *p:source_polygons){
		tree.Insert(p->getMBB()->low, p->getMBB()->high, p);
		p->rasterization(10);
	}
	logt("building R-Tree with %ld nodes", start, source_polygons.size());
	query_context ctx(global_ctx);
	int i=0;
	for(MyPolygon *p:source_polygons){
		struct timeval start = get_cur_time();
		ctx.target = (void *)p;
		Pixel *px = p->getMBB();
		tree.Search(px->low, px->high, MySearchCallback, (void *)&ctx);

		//p->print(false);
		if(++i%1000 == 0){
			logt("queried %d",start,i);
		}
	}

	vector<int> only_children;
	vector<int> only_parent;
	vector<int> both;
	for(MyPolygon *p:source_polygons){
		int id = p->getid();
		bool hasc = has_children.find(id)!=has_children.end();
		bool hasp = has_parent.find(id)!=has_parent.end();
		if(hasc&&!hasp){
			only_parent.push_back(id);
		}
		if(hasp&&!hasc){
			only_children.push_back(id);
		}
		if(hasp&&hasc){
			both.push_back(id);
		}

	}

//	for(int id:only_children){
//		source_polygons[id]->print(false);
//	}
//
//	for(int id:only_parent){
//		source_polygons[id]->print(false);
//	}
//	for(int id:both){
//		source_polygons[id]->print(false);
//	}

	log("%ld only has children %ld only has parent %ld have both",only_parent.size(),only_children.size(),both.size());
	return 0;
}



