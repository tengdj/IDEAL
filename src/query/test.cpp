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
	global_ctx.max_num_polygons = 50;
	vector<MyPolygon *> source_polygons = MyPolygon::load_binary_file("/home/teng/git/IDEAL/data/source.dat",global_ctx, false);

	timeval start = get_cur_time();
	for(MyPolygon *p:source_polygons){
		tree.Insert(p->getMBB()->low, p->getMBB()->high, p);
		p->rasterization(10);
	}
	logt("building R-Tree with %ld nodes", start, source_polygons.size());
	query_context ctx(global_ctx);
	//ctx.geography = true;
	int p1 = 28;
	int p2 = 40;
//	source_polygons[0]->get_rastor()->print();
//	source_polygons[0]->print();
//	log("%d %d %d",source_polygons[0]->get_num_pixels(BORDER),source_polygons[0]->get_num_pixels(IN),source_polygons[0]->get_num_pixels(OUT));
//	return 0;
	int ts = 0;
	int tp = 0;
	int tv = 0;
	int tm = 0;
	int count = 0;
	double right_tm = 0;
	double right_rt = 0;
	double wrong_tm = 0;
	double wrong_rt = 0;
	for(p1=0;p1<source_polygons.size();p1++)
	for(p2=0;p2<source_polygons.size();p2++)
	{
		MyPolygon *target = source_polygons[p1];
		MyPolygon *p = source_polygons[p2];

		if(p->getid()<=target->getid()){
			continue;
		}
		if(p->getid()!=p2){
			//continue;
		}

//		if(p->getMBB()->distance(*target->getMBB(), ctx.geography)>1000){
//			continue;
//		}

		start = get_cur_time();
		double dist = target->distance(p, &ctx);
		double t1 = get_time_elapsed(start, true);
		dist = p->distance(target, &ctx);
		double t2 = get_time_elapsed(start, true);
//		log("%ld %ld %ld %ld",target->get_num_pixels(),target->get_num_pixels(BORDER),
//				p->get_num_pixels(),p->get_num_pixels(BORDER));
		log("id: %d-%d:\ttime: %.2f\t%.2f\tportion: %.3f\t%.3f\t#vertices: %ld\t%ld\tcompare: %d %d %d %d %d",
											target->getid(),p->getid(),
											t1,t2,
											target->get_pixel_portion(IN), p->get_pixel_portion(IN),
											target->get_num_vertices(),p->get_num_vertices(),
											t1>t2,
											target->get_pixel_portion(IN) >p->get_pixel_portion(IN),
											target->get_num_vertices()>p->get_num_vertices(),
											target->get_rastor()->get_step(false)>p->get_rastor()->get_step(false),
											target->getMBB()->area()>p->getMBB()->area()
											);
		count++;
		tp += (t1>t2 == target->get_pixel_portion(IN) >p->get_pixel_portion(IN));
		tv += (t1>t2 == target->get_num_vertices()>p->get_num_vertices());
		ts += (t1>t2 == target->get_rastor()->get_step(false)>p->get_rastor()->get_step(false));
		tm += (t1>t2 == target->getMBB()->area()>p->getMBB()->area());
		if(t1>t2 == target->get_rastor()->get_step(false)>p->get_rastor()->get_step(false)){
			right_tm += max(t1,t2);
			right_rt += max(t1,t2)/min(t1,t2);
		}else{
			wrong_tm += max(t1,t2);
			wrong_rt += max(t1,t2)/min(t1,t2);
		}
		//source_polygons[p1]->print(false);
		//source_polygons[p2]->print(false);
	}
	log("%d %d %d %d %d",count,tp,tv,ts,tm);
	log("%f %f %f %f", right_tm/count,right_rt/count,wrong_tm/count,wrong_rt/count);


//	for(MyPolygon *p:source_polygons){
//		log("%d %d", p->getid(), p->get_num_vertices());
//	}
	return 0;
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



