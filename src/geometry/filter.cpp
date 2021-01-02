/*
 * filter.cpp
 *
 *  Created on: Dec 24, 2020
 *      Author: teng
 */

#include "MyPolygon.h"
#include "../triangulate/poly2tri.h"
using namespace p2t;


void MyPolygon::triangulate(){
	assert(!boundary->clockwise());

	if(!this->valid_for_triangulate){
		return;
	}
	if(rtree){
		delete rtree;
		rtree = NULL;
	}
	if(cdt){
		delete cdt;
		cdt = NULL;
	}
	vector<Point *> polyline = boundary->pack_to_polyline();
	assert(polyline.size()>0);
	cdt = new CDT(polyline);
	cdt->Triangulate();
}

RTree<Triangle *, double, 2, double> * MyPolygon::build_rtree(){
	if(!this->valid_for_triangulate){
		return NULL;
	}
	if(!cdt){
		triangulate();
	}
	assert(cdt);
	if(rtree){
		delete rtree;
		rtree = NULL;
	}
	if(rtree_pixel){
		delete rtree_pixel;
		rtree_pixel = NULL;
	}

	rtree = new RTree<Triangle *, double, 2, double>();
	for(Triangle *tri:cdt->GetTriangles()){
		Pixel pix;
		for(int i=0;i<3;i++){
			pix.update(*tri->point(i));
		}
		rtree->Insert(pix.low, pix.high, tri);
	}
	mbr = getMBB();
	rtree_pixel = new Pixel(mbr);
	rtree->construct_pixel(rtree_pixel);

	return rtree;

}


Pixel *MyPolygon::generateMER(int cx, int cy){
	assert(partitions[cx][cy].status==IN);
	Pixel *curmer = new Pixel();
	int shift[4] = {0,0,0,0};
	bool limit[4] = {false,false,false,false};
	int index = 0;
	while(!limit[0]||!limit[1]||!limit[2]||!limit[3]){
		//left
		if(!limit[0]){
			shift[0]++;
			if(cx-shift[0]<0){
				limit[0] = true;
				shift[0]--;
			}else{
				for(int i=cy-shift[1];i<=cy+shift[3];i++){
					if(partitions[cx-shift[0]][i].status!=IN){
						limit[0] = true;
						shift[0]--;
						break;
					}
				}
			}
		}
		//bottom
		if(!limit[1]){
			shift[1]++;
			if(cy-shift[1]<0){
				limit[1] = true;
				shift[1]--;
			}else{
				for(int i=cx-shift[0];i<=cx+shift[2];i++){
					if(partitions[i][cy-shift[1]].status!=IN){
						limit[1] = true;
						shift[1]--;
						break;
					}
				}
			}
		}
		//right
		if(!limit[2]){
			shift[2]++;
			if(cx+shift[2]>=partitions.size()){
				limit[2] = true;
				shift[2]--;
			}else{
				for(int i=cy-shift[1];i<=cy+shift[3];i++){
					if(partitions[cx+shift[2]][i].status!=IN){
						limit[2] = true;
						shift[2]--;
						break;
					}
				}
			}
		}
		//top
		if(!limit[3]){
			shift[3]++;
			if(cy+shift[3]>=partitions[0].size()){
				limit[3] = true;
				shift[3]--;
			}else{
				for(int i=cx-shift[0];i<=cx+shift[2];i++){
					if(partitions[i][cy+shift[3]].status!=IN){
						limit[3] = true;
						shift[3]--;
						break;
					}
				}
			}
		}
	}

	curmer->low[0] = partitions[cx-shift[0]][cy-shift[1]].low[0];
	curmer->low[1] = partitions[cx-shift[0]][cy-shift[1]].low[1];
	curmer->high[0] = partitions[cx+shift[2]][cy+shift[3]].high[0];
	curmer->high[1] = partitions[cx+shift[2]][cy+shift[3]].high[1];

	return curmer;
}

Pixel *MyPolygon::getMER(query_context *ctx){

	if(mer){
		return mer;
	}
	if(!is_grid_partitioned()){
		partition(ctx->vpr);
	}

	assert(this->is_grid_partitioned());
	vector<int> interiors;
	for(int i=0;i<partitions.size();i++){
		for(int j=0;j<partitions[0].size();j++){
			if(partitions[i][j].status==IN){
				interiors.push_back(i*partitions[0].size()+j);
			}
		}
	}
	if(interiors.size()==0){
		return NULL;
	}
	int loops = ctx->mer_sample_round;
	int maxcount = 0;
	Pixel *max_mer = NULL;
	while(loops-->0){
		int sample = get_rand_number(interiors.size())-1;
		int cx = interiors[sample]/partitions[0].size();
		int cy = interiors[sample]%partitions[0].size();
		Pixel *curmer = generateMER(cx,cy);
		if(max_mer){
			if(max_mer->area()<curmer->area()){
				delete max_mer;
				max_mer = curmer;
			}else{
				delete curmer;
			}
		}else{

			max_mer = curmer;
		}
	}
	interiors.size();
	mer = max_mer;
	return mer;

}

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise
inline int orientation(Point p, Point q, Point r){
	double val = (q.y - p.y) * (r.x - q.x) -
			  (q.x - p.x) * (r.y - q.y);
	if (val == 0.0) return 0;  // colinear
	return (val > 0)? 1: 2; // clock or counterclock wise
}

inline int orientation(double px, double py, double qx, double qy, double rx, double ry){
	double val = (qy - py) * (rx - qx) -
			  (qx - px) * (ry - qy);
	if (val == 0.0) return 0;  // colinear
	return (val > 0)? 1: 2; // clock or counterclock wise
}

// Prints convex hull of a set of n points.
VertexSequence *convexHull(VertexSequence *boundary)
{

	int n = boundary->num_vertices-1;
	assert(n>0);
    // There must be at least 3 points
    if (n < 3){
    	return NULL;
    }
    // Initialize Result
    vector<int> hull;

    // Find the leftmost point
    int p = 0;
    for (int i = 1; i < n; i++)
        if (boundary->x[i] < boundary->x[p])
        	p = i;

    // Start from leftmost point, keep moving counterclockwise
    // until reach the start point again.  This loop runs O(h)
    // times where h is number of points in result or output.
    int idx = 0;
    bool complete = false;
    do
    {
        // Add current point to result
        hull.push_back(p);

        // Search for a point 'q' such that orientation(p, x,
        // q) is counterclockwise for all points 'x'. The idea
        // is to keep track of last visited most counterclock-
        // wise point in q. If any point 'i' is more counterclock-
        // wise than q, then update q.
        int q = (p+1)%n;
        for (int i = 0; i < n; i++)
        {
           // If i is more counterclockwise than current q, then
           // update q
           if (q!=i && orientation(boundary->x[p], boundary->y[p],boundary->x[i], boundary->y[i],boundary->x[q],boundary->y[q]) == 2)
               q = i;
        }

        // Now q is the most counterclockwise with respect to p
        // Set p as q for next iteration, so that q is added to
        // result 'hull'
        p = q;
////
//        cout<<idx++<<"\t"<<p<<"\t"<<left<<"\t"<<boundary->num_vertices;
//        printf("-%f %f-%f %f\n",boundary->x[p],boundary->y[p],boundary->x[left],boundary->y[left]);
//        if(idx++>boundary->num_vertices){
//        	boundary->print();
//        	exit(0);
//        }
        //
		for (int ih=0;ih<hull.size()&&!complete;ih++){
			complete = (p==hull[ih]);
		}
    } while (!complete);  // While we don't come to first point

   VertexSequence *ch = new VertexSequence(hull.size()+1);
   for (int i=0;i<hull.size();i++){
	   ch->x[i] = boundary->x[hull[i]];
	   ch->y[i] = boundary->y[hull[i]];
   }
   ch->x[hull.size()] = ch->x[0];
   ch->y[hull.size()] = ch->y[0];
   hull.clear();
   return ch;
}

VertexSequence *MyPolygon::get_convex_hull(){
	if(convex_hull==NULL){
		convex_hull = convexHull(boundary);
	}
	return convex_hull;
}


void *convex_hull_unit(void *args){

	query_context *ctx = (query_context *)args;
	query_context *gctx = ctx->global_ctx;
	//log("thread %d is started",ctx->thread_id);
	int local_count = 0;
	while(ctx->next_batch(100)){
		for(int i=ctx->index;i<ctx->index_end;i++){
			log("processing\t%d %d",i,gctx->source_polygons[i]->getid());
			VertexSequence *ch = gctx->source_polygons[i]->get_convex_hull();
			log("processed\t%d %d",i,gctx->source_polygons[i]->getid());
			ctx->report_progress();
		}
	}
	ctx->merge_global();
	return NULL;
}

void process_convex_hull(query_context *gctx){

	// must be non-empty
	assert(gctx->source_polygons.size()>0);

	gctx->index = 0;
	size_t former = gctx->target_num;
	gctx->target_num = gctx->source_polygons.size();

	struct timeval start = get_cur_time();
	pthread_t threads[gctx->num_threads];
	query_context ctx[gctx->num_threads];
	for(int i=0;i<gctx->num_threads;i++){
		ctx[i] = *gctx;
		ctx[i].thread_id = i;
		ctx[i].global_ctx = gctx;
	}

	for(int i=0;i<gctx->num_threads;i++){
		pthread_create(&threads[i], NULL, convex_hull_unit, (void *)&ctx[i]);
	}

	for(int i = 0; i < gctx->num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}

	//collect convex hull status
	size_t num_vertexes = 0;
	size_t data_size = 0;
	size_t ch_size = 0;
	for(MyPolygon *poly:gctx->source_polygons){
		num_vertexes += poly->get_convex_hull()->num_vertices;
		data_size += poly->get_data_size();
		ch_size += poly->get_convex_hull()->num_vertices*16;
	}

	logt("get convex hull for %d polygons with %ld average vertices %f overhead", start,
			gctx->source_polygons.size(),
			num_vertexes/gctx->source_polygons.size(),
			ch_size*100.0/data_size);
	gctx->index = 0;
	gctx->query_count = 0;
	gctx->target_num = former;
}

void *mer_unit(void *args){

	query_context *ctx = (query_context *)args;
	query_context *gctx = ctx->global_ctx;
	//log("thread %d is started",ctx->thread_id);
	int local_count = 0;
	while(ctx->next_batch(100)){
		for(int i=ctx->index;i<ctx->index_end;i++){
			gctx->source_polygons[i]->getMER(ctx);
			ctx->report_progress();
		}
	}
	ctx->merge_global();
	return NULL;
}

void process_mer(query_context *gctx){

	// must be non-empty
	assert(gctx->source_polygons.size()>0);

	gctx->index = 0;
	size_t former = gctx->target_num;
	gctx->target_num = gctx->source_polygons.size();

	struct timeval start = get_cur_time();
	pthread_t threads[gctx->num_threads];
	query_context ctx[gctx->num_threads];
	for(int i=0;i<gctx->num_threads;i++){
		ctx[i] = *gctx;
		ctx[i].thread_id = i;
		ctx[i].global_ctx = gctx;
	}

	for(int i=0;i<gctx->num_threads;i++){
		pthread_create(&threads[i], NULL, mer_unit, (void *)&ctx[i]);
	}

	for(int i = 0; i < gctx->num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}

	//collect convex hull status
	double mbr_are = 0;
	double mer_are = 0;
	size_t data_size = 0;
	size_t mer_size = 0;
	for(MyPolygon *poly:gctx->source_polygons){
		data_size += poly->get_data_size();
		mer_size += 4*8;
		mbr_are += poly->getMBB()->area();
		if(poly->getMER(gctx)){
			mer_are += poly->getMER(gctx)->area();
		}
	}

	logt("get mer for %d polygons with %f mer/mbr with %f overhead", start,
			gctx->source_polygons.size(),
			mer_are/mbr_are,
			mer_size*100.0/data_size);
	gctx->index = 0;
	gctx->query_count = 0;
	gctx->target_num = former;
}


void *triangulate_unit(void *args){

	query_context *ctx = (query_context *)args;
	query_context *gctx = ctx->global_ctx;
	//log("thread %d is started",ctx->thread_id);
	int local_count = 0;
	while(ctx->next_batch(100)){
		for(int i=ctx->index;i<ctx->index_end;i++){
			vector<Point*> polyline;
			MyPolygon *polygon = gctx->source_polygons[i];
			if(!polygon->valid_for_triangulate){
				continue;
			}
			polygon->triangulate();
			ctx->report_progress();
		}
	}
	ctx->merge_global();
	return NULL;
}

void process_triangulate(query_context *gctx){

	// must be non-empty
	assert(gctx->source_polygons.size()>0);

	gctx->index = 0;
	size_t former = gctx->target_num;
	gctx->target_num = gctx->source_polygons.size();

	struct timeval start = get_cur_time();
	pthread_t threads[gctx->num_threads];
	query_context ctx[gctx->num_threads];
	for(int i=0;i<gctx->num_threads;i++){
		ctx[i] = *gctx;
		ctx[i].thread_id = i;
		ctx[i].global_ctx = gctx;
	}

	for(int i=0;i<gctx->num_threads;i++){
		pthread_create(&threads[i], NULL, triangulate_unit, (void *)&ctx[i]);
	}

	for(int i = 0; i < gctx->num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}

	//collect convex hull status
	size_t data_size = 0;
	size_t triangle_size = 0;
	for(MyPolygon *poly:gctx->source_polygons){
		data_size += poly->get_data_size();
		triangle_size += poly->get_triangle_size();
	}

	logt("triangulated %d polygons with %f overhead", start,
			gctx->source_polygons.size(),
			triangle_size*100.0/data_size);
	gctx->index = 0;
	gctx->query_count = 0;
	gctx->target_num = former;
}



void *internal_rtree_unit(void *args){

	query_context *ctx = (query_context *)args;
	query_context *gctx = ctx->global_ctx;
	//log("thread %d is started",ctx->thread_id);
	int local_count = 0;
	while(ctx->next_batch(100)){
		for(int i=ctx->index;i<ctx->index_end;i++){
			vector<Point*> polyline;
			MyPolygon *polygon = gctx->source_polygons[i];
			if(!polygon->valid_for_triangulate){
				continue;
			}
			polygon->build_rtree();
			ctx->report_progress();
		}
	}
	ctx->merge_global();
	return NULL;
}

void process_internal_rtree(query_context *gctx){

	// must be non-empty
	assert(gctx->source_polygons.size()>0);

	gctx->index = 0;
	size_t former = gctx->target_num;
	gctx->target_num = gctx->source_polygons.size();

	struct timeval start = get_cur_time();
	pthread_t threads[gctx->num_threads];
	query_context ctx[gctx->num_threads];
	for(int i=0;i<gctx->num_threads;i++){
		ctx[i] = *gctx;
		ctx[i].thread_id = i;
		ctx[i].global_ctx = gctx;
	}

	for(int i=0;i<gctx->num_threads;i++){
		pthread_create(&threads[i], NULL, internal_rtree_unit, (void *)&ctx[i]);
	}

	for(int i = 0; i < gctx->num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}

	//collect convex hull status
	size_t data_size = 0;
	size_t rtree_size = 0;
	for(MyPolygon *poly:gctx->source_polygons){
		data_size += poly->get_data_size();
		rtree_size += poly->get_rtree_size();
	}

	logt("RTree of %d polygons with %f overhead", start,
			gctx->source_polygons.size(),
			rtree_size*100.0/data_size);

	gctx->index = 0;
	gctx->query_count = 0;
	gctx->target_num = former;
}

void preprocess(query_context *gctx){
	if(gctx->use_grid||gctx->use_qtree){
		process_partition(gctx);
	}

	if(gctx->use_convex_hull){
		process_convex_hull(gctx);
	}

	if(gctx->use_mer){
		process_mer(gctx);
	}


	if(gctx->use_triangulate){
		if(gctx->valid_path.size()>0){
			 ifstream is(gctx->valid_path);
			 int num = 0;
			 while(is>>num){
				 assert(num<gctx->source_polygons.size());
				 gctx->source_polygons[num]->valid_for_triangulate = true;
			 }
			 is.close();
		}
		process_triangulate(gctx);
		process_internal_rtree(gctx);
	}
}

