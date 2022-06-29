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

VertexSequence *MyPolygon::get_convex_hull(){
	if(convex_hull==NULL){
		convex_hull = boundary->convexHull();
	}
	return convex_hull;
}
