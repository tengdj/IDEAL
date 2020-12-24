/*
 * filter.cpp
 *
 *  Created on: Dec 24, 2020
 *      Author: teng
 */

#include "MyPolygon.h"


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


// A utility function to find next to top in a stack
Point nextToTop(stack<Point> &S){
	Point p = S.top();
	S.pop();
	Point res = S.top();
	S.push(p);
	return res;
}

// A utility function to swap two points
inline void swap(Point &p1, Point &p2){
	Point temp = p1;
	p1 = p2;
	p2 = temp;
}

// A utility function to return square of distance
// between p1 and p2
inline int distSq(Point p1, Point p2){
	return (p1.x - p2.x)*(p1.x - p2.x) +
		  (p1.y - p2.y)*(p1.y - p2.y);
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
//
//
//map<int, Point> p0_map;
//Point p0;
//
//
//int compare(const void *vp1, const void *vp2)
//{
//   //Point p0 = p0_map[syscall(__NR_gettid)];
//
//   Point *p1 = (Point *)vp1;
//   Point *p2 = (Point *)vp2;
//
//   // Find orientation
//   int o = orientation(p0, *p1, *p2);
//   if (o == 0)
//	 return (distSq(p0, *p2) >= distSq(p0, *p1))? -1 : 1;
//
//   return (o == 2)? -1: 1;
//}
//
//// A function used by library function qsort() to sort an array of
//// points with respect to the first point
//
//int iii = 0;
//int mycompare(Point *p0, Point *p1, Point *p2)
//{
//   // Find orientation
//   int o = orientation(*p0, *p1, *p2);
//   return (o == 2)? -1: 1;
//
//   if (o == 0){
//	   cout<<"coliner "<<iii++<<endl;
//		 return (distSq(*p0, *p2) >= distSq(*p0, *p1))? -1 : 1;
//   }
//
//   return (o == 2)? -1: 1;
//}
//
//
//// Prints convex hull of a set of n points.
//VertexSequence *convexHull2(VertexSequence *boundary){
//	int tid = (int)syscall(__NR_gettid);
//	Point *points = (Point *)malloc(sizeof(Point)*boundary->num_vertices);
//	int n = boundary->num_vertices;
//	for(int i=0;i<boundary->num_vertices;i++){
//		points[i].x = boundary->x[i];
//		points[i].y = boundary->y[i];
//	}
//   // Find the bottommost point
//   int ymin = points[0].y, min = 0;
//   for (int i = 1; i < n; i++)
//   {
//	 int y = points[i].y;
//
//	 // Pick the bottom-most or chose the left
//	 // most point in case of tie
//	 if ((y < ymin) || (ymin == y &&
//		 points[i].x < points[min].x))
//		ymin = points[i].y, min = i;
//   }
//
//   // Place the bottom-most point at first position
//   swap(points[0], points[min]);
//
//   p0_map[tid] = points[0];
//   p0 = points[0];
//
//   // Sort n-1 points with respect to the first point.
//   // A point p1 comes before p2 in sorted output if p2
//   // has larger polar angle (in counterclockwise
//   // direction) than p1
//   qsort(&points[1], n-1, sizeof(Point), compare);
////   for(int i=1;i<boundary->num_vertices-1;i++){
////	   for(int j=i+1;j<boundary->num_vertices;j++){
////		   if(mycompare(&points[0],&points[i],&points[j])>0){
////			   Point tmpp = points[i];
////			   points[i] = points[j];
////			   points[j] = tmpp;
////		   }
////	   }
////   }
//
//   // If two or more points make same angle with p0,
//   // Remove all but the one that is farthest from p0
//   // Remember that, in above sorting, our criteria was
//   // to keep the farthest point at the end when more than
//   // one points have same angle.
//   int m = 1; // Initialize size of modified array
//
//   for (int i=1; i<n; i++)
//   {
//	   // Keep removing i while angle of i and i+1 is same
//	   // with respect to p0
//	   while (i < n-1 && orientation(points[0], points[i], points[i+1]) == 0){
//		   i++;
//	   }
//	   points[m] = points[i];
//	   m++;  // Update size of modified array
//   }
//
//   // If modified array of points has less than 3 points,
//   // convex hull is not possible
//   if (m < 3){
//	   return NULL;
//   }
//
//   // Create an empty stack and push first three points
//   // to it.
//   stack<Point> S;
//   S.push(points[0]);
//   S.push(points[1]);
//   S.push(points[2]);
//
//   // Process remaining n-3 points
//   for (int i = 3; i < m; i++)
//   {
//	  // Keep removing top while the angle formed by
//	  // points next-to-top, top, and points[i] makes
//	  // a non-left turn
//	  while (orientation(nextToTop(S), S.top(), points[i]) != 2)
//		 S.pop();
//	  S.push(points[i]);
//   }
//
//   VertexSequence *ch = new VertexSequence();
//   ch->x = new double[S.size()+1];
//   ch->y = new double[S.size()+1];
//
//   // Now stack has the output points
//   while (!S.empty())
//   {
//	   Point p = S.top();
//	   ch->x[ch->num_vertices] = p.x;
//	   ch->y[ch->num_vertices++] = p.y;
//	   S.pop();
//   }
//   ch->x[ch->num_vertices] = ch->x[0];
//   ch->y[ch->num_vertices++] = ch->y[0];
//
//   free(points);
//   return ch;
//}



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
