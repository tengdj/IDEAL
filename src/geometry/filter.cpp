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
	for(Triangle *tri:triangles){
		delete tri;
	}
	triangles.clear();

	vector<Vertex *> polyline = boundary->pack_to_polyline();
	assert(polyline.size()>0);
	CDT *cdt = new CDT(polyline);
	cdt->Triangulate();
	vector<Triangle *> tri = cdt->GetTriangles();
	triangles.resize(tri.size());
	for(int i=0;i<tri.size();i++){
		triangles[i] = new Triangle(*tri[i]->point(0),*tri[i]->point(1),*tri[i]->point(2));
	}

	for(Vertex *v:polyline){
		delete v;
	}
	polyline.clear();
	delete cdt;
}

void MyPolygon::build_rtree(){
	if(!this->valid_for_triangulate){
		return;
	}

	triangulate();

	if(rtree){
		delete rtree;
		rtree = NULL;
	}

	RTree<Triangle *, double, 2, double> *rtree_tmp = new RTree<Triangle *, double, 2, double>();
	for(Triangle *tri:triangles){
		Pixel pix;
		for(int i=0;i<3;i++){
			pix.update(*tri->point(i));
		}
		rtree_tmp->Insert(pix.low, pix.high, tri);
	}
	mbr = getMBB();
	rtree = new RTNode();
	rtree->low[0] = mbr->low[0];
	rtree->low[1] = mbr->low[1];
	rtree->high[0] = mbr->high[0];
	rtree->high[1] = mbr->high[1];

	rtree_tmp->construct_pixel(rtree);
	delete rtree_tmp;
}

Pixel *MyPolygon::getMBB(){
	if(this->mbr){
		return mbr;
	}
	mbr = boundary->getMBR();
	return mbr;
}

Pixel *MyPolygon::getMER(query_context *ctx){

	if(mer){
		return mer;
	}

	MyRaster *ras = new MyRaster(boundary,ctx->vpr);
	vector<Pixel *> interiors = ras->get_pixels(IN);
	if(interiors.size()==0){
		return NULL;
	}
	int loops = ctx->mer_sample_round;
	int maxcount = 0;
	Pixel *max_mer = NULL;
	while(loops-->0){
		int sample = get_rand_number(interiors.size())-1;
		Pixel *curmer = raster->extractMER(interiors[sample]);
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
	interiors.clear();
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
    // There must be at least 3 points
    if (n < 3){
    	return NULL;
    }
    // Initialize Result
    vector<int> hull;

    // Find the leftmost point
    int p = 0;
    for (int i = 1; i < n; i++)
        if (boundary->p[i].x < boundary->p[p].x)
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
           if (q!=i && orientation(boundary->p[p].x, boundary->p[p].y,boundary->p[i].x, boundary->p[i].y,boundary->p[q].x,boundary->p[q].y) == 2)
               q = i;
        }

        // Now q is the most counterclockwise with respect to p
        // Set p as q for next iteration, so that q is added to
        // result 'hull'
        p = q;
		for (int ih=0;ih<hull.size()&&!complete;ih++){
			complete = (p==hull[ih]);
		}
    } while (!complete);  // While we don't come to first point

   VertexSequence *ch = new VertexSequence(hull.size()+1);
   for (int i=0;i<hull.size();i++){
	   ch->p[i].x = boundary->p[hull[i]].x;
	   ch->p[i].y = boundary->p[hull[i]].y;
   }
   ch->p[hull.size()].x = ch->p[0].x;
   ch->p[hull.size()].y = ch->p[0].y;
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
			gctx->source_polygons[i]->get_convex_hull();
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
		VertexSequence *ch = poly->get_convex_hull();
		if(ch){
			num_vertexes += ch->num_vertices;
			ch_size += ch->num_vertices*16;
		}
		data_size += poly->get_data_size();
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
		process_rasterization(gctx);
	}

	if(gctx->use_convex_hull){
		process_convex_hull(gctx);
	}

	if(gctx->use_mer){
		process_mer(gctx);
	}

	int idx = 0;
	if(gctx->use_triangulate){
		if(gctx->valid_path.size()>0){
			 ifstream is(gctx->valid_path);
			 int num = 0;
			 while(is>>num){
				 assert(num<gctx->source_polygons.size());
				 gctx->source_polygons[num]->valid_for_triangulate = gctx->source_polygons[num]->area()>0&&true;
			 }
			 is.close();
		}
		process_triangulate(gctx);
		process_internal_rtree(gctx);
	}
}




void *partition_unit(void *args){
	query_context *ctx = (query_context *)args;
	query_context *gctx = ctx->global_ctx;
	//log("thread %d is started",ctx->thread_id);
	int local_count = 0;
	while(ctx->next_batch(1)){
		for(int i=ctx->index;i<ctx->index_end;i++){
			struct timeval start = get_cur_time();

			if(ctx->use_grid){
				gctx->source_polygons[i]->rasterization(ctx->vpr);
			}
			if(ctx->use_qtree){
				gctx->source_polygons[i]->partition_qtree(ctx->vpr);
			}
			double latency = get_time_elapsed(start);
			int num_vertices = gctx->source_polygons[i]->get_num_vertices();
			//ctx->report_latency(num_vertices, latency);
			if(latency>10000||num_vertices>200000){
				logt("partition %d vertices (source)",start,num_vertices);
			}
			ctx->report_progress();
		}
	}
	ctx->merge_global();
	return NULL;
}

void *partition_target_unit(void *args){
	query_context *ctx = (query_context *)args;
	query_context *gctx = ctx->global_ctx;
	//log("thread %d is started",ctx->thread_id);
	int local_count = 0;
	while(ctx->next_batch(1)){
		for(int i=ctx->index;i<ctx->index_end;i++){
			struct timeval start = get_cur_time();

			if(ctx->use_grid){
				gctx->target_polygons[i]->rasterization(ctx->vpr);
			}
			if(ctx->use_qtree){
				gctx->target_polygons[i]->partition_qtree(ctx->vpr);
			}
			double latency = get_time_elapsed(start);
			int num_vertices = gctx->target_polygons[i]->get_num_vertices();
			//ctx->report_latency(num_vertices, latency);
//			if(latency>10000||num_vertices>200000){
//				logt("partition %d vertices (target)",start,num_vertices);
//			}
			ctx->report_progress();
		}
	}
	ctx->merge_global();
	return NULL;
}

void process_rasterization(query_context *gctx){

	log("start rasterizing the referred polygons");
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
		pthread_create(&threads[i], NULL, partition_unit, (void *)&ctx[i]);
	}

	for(int i = 0; i < gctx->num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}

	//collect partitioning status
	size_t num_partitions = 0;
	size_t num_crosses = 0;
	size_t num_border_partitions = 0;
	size_t num_edges = 0;
	for(MyPolygon *poly:gctx->source_polygons){
		if(gctx->use_grid){
			num_partitions += poly->get_rastor()->get_num_pixels();
			num_crosses += poly->get_rastor()->get_num_crosses();
			num_border_partitions += poly->get_rastor()->get_num_pixels(BORDER);
			num_edges += poly->get_rastor()->get_num_border_edge();
		}else if(gctx->use_qtree){
			num_partitions += poly->get_qtree()->leaf_count();
			num_border_partitions += poly->get_qtree()->border_leaf_count();
		}
	}
	logt("partitioned %d polygons with (%ld)%ld average pixels %.2f average crosses per pixel %.2f edges per pixel", start,
			gctx->source_polygons.size(),
			num_border_partitions/gctx->source_polygons.size(),
			num_partitions/gctx->source_polygons.size(),
			1.0*num_crosses/num_border_partitions,
			1.0*num_edges/num_border_partitions);

	if(gctx->target_polygons.size()>0){

		gctx->index = 0;
		gctx->target_num = gctx->target_polygons.size();
		gctx->query_count = 0;

		struct timeval start = get_cur_time();
		pthread_t threads[gctx->num_threads];
		query_context ctx[gctx->num_threads];
		for(int i=0;i<gctx->num_threads;i++){
			ctx[i] = *gctx;
			ctx[i].thread_id = i;
			ctx[i].global_ctx = gctx;
		}
		for(int i=0;i<gctx->num_threads;i++){
			pthread_create(&threads[i], NULL, partition_target_unit, (void *)&ctx[i]);
		}

		for(int i = 0; i < gctx->num_threads; i++ ){
			void *status;
			pthread_join(threads[i], &status);
		}

		size_t num_partitions = 0;
		size_t num_crosses = 0;
		size_t num_border_partitions = 0;
		size_t num_edges = 0;
		for(MyPolygon *poly:gctx->target_polygons){
			if(gctx->use_grid){
				num_partitions += poly->get_rastor()->get_num_pixels();
				num_crosses += poly->get_rastor()->get_num_crosses();
				num_border_partitions += poly->get_rastor()->get_num_pixels(BORDER);
				num_edges += poly->get_rastor()->get_num_border_edge();
			}else if(gctx->use_qtree){
				num_partitions += poly->get_qtree()->leaf_count();
				num_border_partitions += poly->get_qtree()->border_leaf_count();
			}
		}
		logt("partitioned %d polygons with (%ld)%ld average pixels %.2f average crosses per pixel %.2f edges per pixel (target)", start,
				gctx->target_polygons.size(),
				num_border_partitions/gctx->target_polygons.size(),
				num_partitions/gctx->target_polygons.size(),
				1.0*num_crosses/num_border_partitions,
				1.0*num_edges/num_border_partitions);
	}


	gctx->index = 0;
	gctx->query_count = 0;
	gctx->target_num = former;
}




