/*
 * filter.cpp
 *
 *  Created on: Dec 24, 2020
 *      Author: teng
 */

#include "MyPolygon.h"
#include "../triangulate/poly2tri.h"
using namespace p2t;


int  *MyPolygon::triangulate(){
	assert(!boundary->clockwise());
	triangles = polygon_triangulate(boundary->num_vertices-1, boundary->x, boundary->y);
	return triangles;
}

void MyPolygon::build_rtree(){
	if(!triangles){
		triangulate();
	}
	assert(triangles);
	if(rtree){
		delete rtree;
	}
	rtree = new RTree<int *, double, 2, double>();

	for(int n=0;n<boundary->num_vertices-2;n++){
		double low[2] = {boundary->x[triangles[3*n]], boundary->y[triangles[3*n]]};
		double high[2] = {boundary->x[triangles[3*n]], boundary->y[triangles[3*n]]};
		for(int i=1;i<3;i++){
			double cur_x = boundary->x[triangles[i+3*n]];
			if(cur_x<low[0]){
				low[0] = cur_x;
			}else if(cur_x>high[0]){
				high[0] = cur_x;
			}
			double cur_y = boundary->y[triangles[i+3*n]];
			if(cur_y<low[1]){
				low[1] = cur_y;
			}else if(cur_y>high[1]){
				high[1] = cur_y;
			}
		}
		rtree->Insert(low, high, triangles+3*n);
	}

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

// Prints convex hull of a set of n points.
VertexSequence *convexHull(VertexSequence *boundary)
{
	Point *points = (Point *)malloc(sizeof(Point)*boundary->num_vertices);
	int n = boundary->num_vertices;
	for(int i=0;i<boundary->num_vertices;i++){
		points[i].x = boundary->x[i];
		points[i].y = boundary->y[i];
	}
    // There must be at least 3 points
    if (n < 3){
    	return NULL;
    }

    // Initialize Result
    vector<Point> hull;

    // Find the leftmost point
    int l = 0;
    for (int i = 1; i < n; i++)
        if (points[i].x < points[l].x)
            l = i;

    // Start from leftmost point, keep moving counterclockwise
    // until reach the start point again.  This loop runs O(h)
    // times where h is number of points in result or output.
    int p = l, q;
    do
    {
        // Add current point to result
        hull.push_back(points[p]);

        // Search for a point 'q' such that orientation(p, x,
        // q) is counterclockwise for all points 'x'. The idea
        // is to keep track of last visited most counterclock-
        // wise point in q. If any point 'i' is more counterclock-
        // wise than q, then update q.
        q = (p+1)%n;
        for (int i = 0; i < n; i++)
        {
           // If i is more counterclockwise than current q, then
           // update q
           if (q!=i && orientation(points[p], points[i], points[q]) == 2)
               q = i;
        }

        // Now q is the most counterclockwise with respect to p
        // Set p as q for next iteration, so that q is added to
        // result 'hull'
        p = q;

    } while (p != l);  // While we don't come to first point

   VertexSequence *ch = new VertexSequence();
   ch->x = new double[hull.size()+1];
   ch->y = new double[hull.size()+1];

   for (int i=0;i<hull.size();i++){
	   ch->x[ch->num_vertices] = hull[i].x;
	   ch->y[ch->num_vertices++] = hull[i].y;
   }
   ch->x[ch->num_vertices] = ch->x[0];
   ch->y[ch->num_vertices++] = ch->y[0];

   free(points);
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
			struct timeval start = get_cur_time();
			VertexSequence *ch = gctx->source_polygons[i]->get_convex_hull();
			double latency = get_time_elapsed(start);
			int num_vertices = gctx->source_polygons[i]->get_num_vertices();
			if(latency>10000){
				logt("get convex hull from polygon with %d vertices",start,num_vertices);
			}
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
	for(MyPolygon *poly:gctx->source_polygons){
		num_vertexes += poly->get_convex_hull()->num_vertices;
	}

	logt("get convex hull for %d polygons with %ld average vertices", start,
			gctx->source_polygons.size(),
			num_vertexes/gctx->source_polygons.size());
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
	double mbr_size = 0;
	double mer_size = 0;
	for(MyPolygon *poly:gctx->source_polygons){
		mbr_size = poly->getMBB()->area();
		if(poly->getMER()){
			mer_size = poly->getMER()->area();
		}
	}

	logt("get mer for %d polygons with %f mer/mbr", start,
			gctx->source_polygons.size(),
			mer_size/mbr_size);
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
			vector<TrPoint*> polyline;
			MyPolygon *polygon = gctx->source_polygons[i];
			if(!polygon->valid_for_triangulate){
				continue;
			}
			for(int i=0;i<polygon->boundary->num_vertices-1;i++){
				polyline.push_back(new TrPoint(polygon->boundary->x[i],polygon->boundary->y[i]));
			}
			struct timeval start = get_cur_time();
			CDT* cdt = new CDT(polyline);
			cdt->Triangulate();
			vector<Triangle*> triangles = cdt->GetTriangles();
			delete cdt;
			polyline.clear();
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
//	double mbr_size = 0;
//	double mer_size = 0;
//	for(MyPolygon *poly:gctx->source_polygons){
//		mbr_size = poly->getMBB()->area();
//		if(poly->getMER()){
//			mer_size = poly->getMER()->area();
//		}
//	}
//
//	logt("triangulated %d polygons with %f mer/mbr", start,
//			gctx->source_polygons.size(),
//			mer_size/mbr_size);
	logt("triangulated %d polygons", start,
			gctx->source_polygons.size());
	gctx->index = 0;
	gctx->query_count = 0;
	gctx->target_num = former;
}

