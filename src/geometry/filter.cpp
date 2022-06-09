/*
 * filter.cpp
 *
 *  Created on: Dec 24, 2020
 *      Author: teng
 */

#include "../include/MyPolygon.h"
#include "../triangulate/poly2tri.h"
using namespace p2t;


void MyPolygon::triangulate(){
	assert(!boundary->clockwise());
	if(triangle_num > 0){
		return;
	}
	vector<Vertex *> polyline = boundary->pack_to_polyline();
	assert(polyline.size() > 0);
	CDT *cdt = new CDT(polyline);
	cdt->Triangulate();
	vector<Triangle *> tri = cdt->GetTriangles();
	triangles = new Point[3*tri.size()];
	for(int i=0;i<tri.size();i++){
		triangles[i*3].x = tri[i]->point(0)->x;
		triangles[i*3].y = tri[i]->point(0)->y;
		triangles[i*3+1].x = tri[i]->point(1)->x;
		triangles[i*3+1].y = tri[i]->point(1)->y;
		triangles[i*3+2].x = tri[i]->point(2)->x;
		triangles[i*3+2].y = tri[i]->point(2)->y;
	}
	triangle_num = tri.size();
	delete cdt;
	for(Vertex *v:polyline){
		delete v;
	}
	polyline.clear();
}

void MyPolygon::build_rtree(){
	triangulate();

	assert(triangles && triangle_num>0);
	if(rtree){
		return;
	}

	RTree<Point *, double, 2, double> *rtree_tmp = new RTree<Point *, double, 2, double>();

	for(int i=0;i<triangle_num;i++){
		box pix;
		Point *ps = triangles+3*i;
		for(int i=0;i<3;i++){
			pix.update(ps[i]);
		}
		rtree_tmp->Insert(pix.low, pix.high, ps);
	}
	mbr = getMBB();
	rtree = new RTNode();
	rtree->low[0] = mbr->low[0];
	rtree->low[1] = mbr->low[1];
	rtree->high[0] = mbr->high[0];
	rtree->high[1] = mbr->high[1];

	rtree_tmp->construct_pixel(rtree);
	assert(rtree->validate());
	delete rtree_tmp;
}

box *MyPolygon::getMBB(){
	if(this->mbr){
		return mbr;
	}
	mbr = boundary->getMBR();
	return mbr;
}

box *MyPolygon::getMER(query_context *ctx){

	if(mer){
		return mer;
	}

	MyRaster *ras = new MyRaster(boundary,ctx->vpr);
	ras->rasterization();
	vector<Pixel *> interiors = ras->get_pixels(IN);
	if(interiors.size()==0){
		return NULL;
	}
	int loops = ctx->mer_sample_round;
	int maxcount = 0;
	box *max_mer = NULL;
	while(loops-->0){
		int sample = get_rand_number(interiors.size())-1;
		box *curmer = ras->extractMER(interiors[sample]);
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
	delete ras;
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
