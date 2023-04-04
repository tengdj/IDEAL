/*
 * Parser.cpp
 *
 *  Created on: May 9, 2020
 *      Author: teng
 */



#include "../index/RTree.h"
#include <queue>
#include <fstream>
#include "../include/MyPolygon.h"
#include <omp.h>
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

	if(poly->contain(target, ctx)){
		has_children.insert(poly->getid());
		has_parent.insert(target->getid());
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


//	printf("ID: %d, Max threads: %d, Num threads: %d\n",omp_get_thread_num(), omp_get_max_threads(), omp_get_num_threads());
//	omp_set_num_threads(5);
//	printf("ID: %d, Max threads: %d, Num threads: %d\n",omp_get_thread_num(), omp_get_max_threads(), omp_get_num_threads());
//
//#pragma omp parallel num_threads(12)
//	{
//		// omp_set_num_threads(6);	// Do not call it in parallel region
//		printf("ID: %d, Max threads: %d, Num threads: %d\n",omp_get_thread_num(), omp_get_max_threads(), omp_get_num_threads());
//	}
//
//	printf("ID: %d, Max threads: %d, Num threads: %d\n",omp_get_thread_num(), omp_get_max_threads(), omp_get_num_threads());
//
//	omp_set_num_threads(6);
//	printf("ID: %d, Max threads: %d, Num threads: %d\n",omp_get_thread_num(), omp_get_max_threads(), omp_get_num_threads());
	query_context cctx;
	cctx.query_type = QueryType::distance;
	cctx.geography = false;


	vector<MyPolygon *> polys = load_binary_file(argv[1], cctx);
	MyPolygon *poly = polys[0];
	//poly->print();
	//Point p(-0.4, 0);
	Point p(0, 0);
	poly->rasterization(10);
	double dist2 = poly->distance(p, &cctx);
	cctx.print_stats();
	double dist1 = poly->boundary->distance(p, false);

	cout<<dist1<<" "<<dist2<<endl;

	//poly->get_rastor()->print();


	return 0;
	query_context global_ctx;
	//global_ctx.max_num_polygons = 200;
	global_ctx.use_grid = true;
	global_ctx.vpr = 10;
	global_ctx.geography = true;
	timeval start = get_cur_time();

	//global_ctx.max_num_polygons = wrong+1;
	vector<MyPolygon *> source_polygons = load_binary_file("/home/teng/git/IDEAL/src/hild_valid.dat",global_ctx);
	global_ctx.source_polygons = source_polygons;
	logt("loading data",start);



	preprocess(&global_ctx);
	logt("preprocess data",start);

	for(MyPolygon *p:source_polygons){
		tree.Insert(p->getMBB()->low, p->getMBB()->high, p);
	}
	logt("building R-Tree with %ld nodes", start, source_polygons.size());
	query_context ctx(global_ctx);

	int ts = 0;
	int tp = 0;
	int tv = 0;
	int tm = 0;
	int count = 0;
	double right_tm = 0;
	double right_rt = 0;
	double wrong_tm = 0;
	double wrong_rt = 0;
	int right_ct = 0;
	int wrong_ct = 0;

	double min_tm = 0;
	double max_tm = 0;

	double choose_big_tm = 0;
	double choose_small_tm = 0;

	for(int p1=0;p1<source_polygons.size()-1;p1++){
		for(int p2=p1+1;p2<source_polygons.size();p2++){
			MyPolygon *target = source_polygons[p1];
			MyPolygon *source = source_polygons[p2];

			start = get_cur_time();
			double dist = target->distance(source, &ctx);
			double t1 = get_time_elapsed(start, true);
			dist = source->distance(target, &ctx);
			double t2 = get_time_elapsed(start, true);
	//		log("%ld %ld %ld %ld",target->get_num_pixels(),target->get_num_pixels(BORDER),
	//				p->get_num_pixels(),p->get_num_pixels(BORDER));
			if(t1>t2 == target->get_rastor()->get_step(false)<source->get_rastor()->get_step(false))
			log("id: %d-%d:\ttime: %.2f\t%.2f\tportion: %.3f\t%.3f\t#vertices: %ld\t%ld\tcompare: %d %d %d %d %d",
												target->getid(),source->getid(),
												t1,t2,
												target->get_pixel_portion(IN), source->get_pixel_portion(IN),
												target->get_num_vertices(),source->get_num_vertices(),
												t1>t2,
												target->get_pixel_portion(IN) >source->get_pixel_portion(IN),
												target->get_num_vertices()>source->get_num_vertices(),
												target->get_rastor()->get_step(false)>source->get_rastor()->get_step(false),
												target->getMBB()->area()>source->getMBB()->area()
												);
			count++;
	//		tp += (t1>t2 == target->get_pixel_portion(IN) >p->get_pixel_portion(IN));
	//		tv += (t1>t2 == target->get_num_vertices()>p->get_num_vertices());
	//		tm += (t1>t2 == target->getMBB()->area()>p->getMBB()->area());

			if(target->get_rastor()->get_step(false)>source->get_rastor()->get_step(false)){
				choose_big_tm += t1;
				choose_small_tm += t2;
			}else{
				choose_big_tm += t2;
				choose_small_tm += t1;
			}

			bool guess_right = ((t1>t2) == (target->get_rastor()->get_step(false)>source->get_rastor()->get_step(false)));
			if(guess_right){
				right_tm += max(t1,t2);
				right_rt += max(t1,t2)/min(t1,t2);
				right_ct++;
			}else{
				wrong_tm += max(t1,t2);
				wrong_rt += max(t1,t2)/min(t1,t2);
				wrong_ct++;
			}
			min_tm += min(t1, t2);
			max_tm += max(t1, t2);
			//source_polygons[p1]->print(false);
			//source_polygons[p2]->print(false);
		}
	}
	//log("%d %d %d %d %d",count,tp,tv,ts,tm);
	log("%d = %d + %d",count,right_ct, wrong_ct);
	log("%f %f %f %f",right_tm/right_ct,right_rt/right_ct,wrong_tm/wrong_ct,wrong_rt/wrong_ct);

	log("%f %f (%f - %f)", choose_big_tm/count, choose_small_tm/count,min_tm/count,max_tm/count);


	if(false){
		int i=0;
		start = get_cur_time();
		for(MyPolygon *p:source_polygons){
			ctx.target = (void *)p;
			box *px = p->getMBB();
			tree.Search(px->low, px->high, MySearchCallback, (void *)&ctx);

			//p->print(false);
			if(++i%1000 == 0){
				logt("queried %d",start,i);
			}
		}

		vector<int> only_children;
		vector<int> only_parent;
		vector<int> both;
		vector<int> neither;
		for(MyPolygon *p:source_polygons){
			int id = p->getid();
			bool hasc = has_children.find(id)!=has_children.end();
			bool hasp = has_parent.find(id)!=has_parent.end();
			if(hasc&&!hasp){
				only_parent.push_back(id);
			}
			if(hasc&&hasp){
				both.push_back(id);
			}
			if(!hasc&&hasp){
				only_children.push_back(id);
			}
			if(!hasc&&!hasp){
				neither.push_back(id);
			}
		}

		for(int id:only_children){
			source_polygons[id]->print(false);
		}
		for(int id:only_parent){
			source_polygons[id]->print(false);
		}
		for(int id:both){
			source_polygons[id]->print(false);
		}
		for(int id:neither){
			source_polygons[id]->print(false);
		}

		log("%ld only are children %ld only are parent %ld are both %ld are neither",only_children.size(),only_parent.size(),both.size(),neither.size());
		return 0;

	}
	return 0;
}



