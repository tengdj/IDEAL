/*
 * MyPolygon.cpp
 *
 *  Created on: May 9, 2020
 *      Author: teng
 */

#include <string.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <float.h>
#include <sstream>
#include <vector>
#include <thread>
#include <unordered_map>

#include "MyPolygon.h"

using namespace std;

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

VertexSequence *VertexSequence::convexHull()
{

	int n = num_vertices-1;
    // There must be at least 3 points
    if (n < 3){
    	return NULL;
    }
    // Initialize Result
    vector<int> hull;

    // Find the leftmost point
    int lp = 0;
    for (int i = 1; i < n; i++)
        if (p[i].x < p[lp].x)
        	lp = i;

    // Start from leftmost point, keep moving counterclockwise
    // until reach the start point again.  This loop runs O(h)
    // times where h is number of points in result or output.
    int idx = 0;
    bool complete = false;
    do
    {
        // Add current point to result
        hull.push_back(lp);

        // Search for a point 'q' such that orientation(lp, x,
        // q) is counterclockwise for all points 'x'. The idea
        // is to keep track of last visited most counterclock-
        // wise point in q. If any point 'i' is more counterclock-
        // wise than q, then update q.
        int q = (lp+1)%n;
        for (int i = 0; i < n; i++)
        {
           // If i is more counterclockwise than current q, then
           // update q
           if (q!=i && orientation(p[lp].x, p[lp].y,p[i].x, p[i].y,p[q].x,p[q].y) == 2)
               q = i;
        }

        // Now q is the most counterclockwise with respect to lp
        // Set lp as q for next iteration, so that q is added to
        // result 'hull'
        lp = q;
		for (int ih=0;ih<hull.size()&&!complete;ih++){
			complete = (lp==hull[ih]);
		}
    } while (!complete);  // While we don't come to first point

   VertexSequence *ch = new VertexSequence(hull.size()+1);
   for (int i=0;i<hull.size();i++){
	   ch->p[i].x = p[hull[i]].x;
	   ch->p[i].y = p[hull[i]].y;
   }
   ch->p[hull.size()].x = ch->p[0].x;
   ch->p[hull.size()].y = ch->p[0].y;
   hull.clear();
   return ch;
}

void VertexSequence::print(bool complete_ring){
	cout<<"(";
	for(int i=0;i<num_vertices;i++){
		if(i!=0){
			cout<<",";
		}
		printf("%f ",p[i].x);
		printf("%f",p[i].y);
	}
	// the last vertex should be the same as the first one for a complete ring
	if(complete_ring){
		if(p[0].x!=p[num_vertices-1].x||p[0].y!=p[num_vertices-1].y){
			cout<<",";
			printf("%f ",p[0].x);
			printf("%f",p[0].y);
		}
	}
	cout<<")";
}

double VertexSequence::area(){
	double sum = 0;
	for(int i=0;i<num_vertices;i++){
		sum += (p[i].x-p[(i+1)%num_vertices].x)*(p[(i+1)%num_vertices].y+p[i].y);
	}
	return sum/2;
}

bool VertexSequence::clockwise(){
	return area()<0;
}

void VertexSequence::reverse(){
	Point tmp;
	for(int i=0;i<num_vertices/2;i++){
		tmp = p[i];
		p[i] = p[num_vertices-1-i];
		p[num_vertices-1-i] = tmp;
	}
}

bool MyPolygon::validate_vertices(const char *wkt, size_t &offset, size_t &len){
	// read until the left parenthesis
	skip_space(wkt, offset);
	if(wkt[offset++]!='('){
		return false;
	}

	// count the number of vertices
	int cur_offset = offset;
	int num_vertices = 0;
	while(wkt[cur_offset++]!=')'){
		if(wkt[cur_offset]==','){
			num_vertices++;
		}
	}
	num_vertices++;

	// read x/y
	for(int i=0;i<num_vertices;i++){
		read_double(wkt, offset);
		read_double(wkt, offset);
	}

	// read until the right parenthesis
	skip_space(wkt, offset);
	if(wkt[offset++]!=')'){
		return false;
	}
	return true;
}


VertexSequence *MyPolygon::read_vertices(const char *wkt, size_t &offset, bool clockwise){
	// read until the left parenthesis
	skip_space(wkt, offset);
	assert(wkt[offset++]=='(');

	// count the number of vertices
	int cur_offset = offset;
	int num_vertices = 0;
	while(wkt[cur_offset++]!=')'){
		if(wkt[cur_offset]==','){
			num_vertices++;
		}
	}
	num_vertices++;
	VertexSequence *vs = new VertexSequence(num_vertices);

	// read x/y
	for(int i=0;i<num_vertices;i++){
		vs->p[i].x = read_double(wkt, offset);
		vs->p[i].y = read_double(wkt, offset);
	}
	if(clockwise){
		if(!vs->clockwise()){
			vs->reverse();
		}
	}else{
		if(vs->clockwise()){
			vs->reverse();
		}
	}

	// read until the right parenthesis
	skip_space(wkt, offset);
	assert(wkt[offset++]==')');
	return vs;
}

MyPolygon *MyPolygon::clone(){
	MyPolygon *polygon = new MyPolygon();
	polygon->boundary = boundary->clone();
	for(VertexSequence *t:holes){
		polygon->holes.push_back(t->clone());
	}
	polygon->getMBB();
	return polygon;
}

bool MyPolygon::validate_polygon(const char *wkt, size_t &offset, size_t &len){
	skip_space(wkt, offset);
	// left parentheses for the entire polygon
	if(wkt[offset++]!='('){
		return false;
	}

	// read the vertices of the boundary polygon
	// the vertex must rotation in clockwise
	if(!validate_vertices(wkt, offset, len)){
		return false;
	}
	skip_space(wkt, offset);
	//polygons as the holes of the boundary polygon
	while(wkt[offset]==','){
		offset++;
		if(!validate_vertices(wkt, offset, len)){
			return false;
		}
		skip_space(wkt, offset);
	}
	if(wkt[offset++]!=')'){
		return false;
	}
	return true;
}

MyPolygon *MyPolygon::read_polygon(const char *wkt, size_t &offset){
	MyPolygon *polygon = new MyPolygon();
	skip_space(wkt, offset);
	// left parentheses for the entire polygon
	assert(wkt[offset++]=='(');

	// read the vertices of the boundary polygon
	// the vertex must rotation in clockwise
	polygon->boundary = read_vertices(wkt, offset,false);
	if(polygon->boundary->clockwise()){
		polygon->boundary->reverse();
	}
	polygon->boundary->fix();
	skip_space(wkt, offset);
	//polygons as the holes of the boundary polygon
	while(wkt[offset]==','){
		offset++;
		VertexSequence *vc = read_vertices(wkt, offset,true);
		if(!vc->clockwise()){
			vc->reverse();
		}
		vc->fix();
		polygon->holes.push_back(vc);

		skip_space(wkt, offset);
	}
	assert(wkt[offset++]==')');
	polygon->getMBB();
	return polygon;
}

MyPolygon::~MyPolygon(){
	clear();
}
void MyPolygon::clear(){
	if(boundary){
		delete boundary;
	}
	for(VertexSequence *p:holes){
		if(p){
			delete p;
		}
	}
	holes.clear();
	if(mbr){
		delete mbr;
	}
	if(mer){
		delete mer;
	}
	if(raster){
		delete raster;
	}
	if(qtree){
		delete qtree;
	}
	if(rtree){
		delete rtree;
	}
	if(triangles){
		delete []triangles;
	}
}

PolygonMeta MyPolygon::get_meta(){
	PolygonMeta pmeta;
	pmeta.size = get_data_size();
	pmeta.num_vertices = get_num_vertices();
	pmeta.mbr = *getMBB();
	return pmeta;
}

size_t MyPolygon::get_data_size(){
	size_t ds = 0;
	ds += sizeof(size_t);
	// for boundary
	ds += boundary->get_data_size();
	for(VertexSequence *vs:holes){
		ds += vs->get_data_size();
	}
	return ds;
}

/*
 * |num_holes|boundary|holes boundaries|
 *
 * */
size_t MyPolygon::encode(char *target){
	size_t encoded = 0;
	((size_t *)target)[0] = holes.size();
	encoded += sizeof(size_t); //saved one size_t for number of holes
	encoded += boundary->encode(target+encoded);
	for(VertexSequence *vs:holes){
		encoded += vs->encode(target+encoded);
	}
	return encoded;
}
size_t MyPolygon::decode(char *source){
	size_t decoded = 0;
	assert(!boundary);
	boundary = new VertexSequence();
	size_t num_holes = ((size_t *)source)[0];
	decoded += sizeof(size_t);
	decoded += boundary->decode(source+decoded);
	for(size_t i=0;i<num_holes;i++){
		VertexSequence *vs = new VertexSequence();
		decoded += vs->decode(source+decoded);
	}
	return decoded;
}

MyPolygon *MyPolygon::gen_box(double min_x,double min_y,double max_x,double max_y){
	MyPolygon *mbr = new MyPolygon();
	mbr->boundary = new VertexSequence(5);
	mbr->boundary->p[0].x = min_x;
	mbr->boundary->p[0].y = min_y;
	mbr->boundary->p[1].x = max_x;
	mbr->boundary->p[1].y = min_y;
	mbr->boundary->p[2].x = max_x;
	mbr->boundary->p[2].y = max_y;
	mbr->boundary->p[3].x = min_x;
	mbr->boundary->p[3].y = max_y;
	mbr->boundary->p[4].x = min_x;
	mbr->boundary->p[4].y = min_y;
	return mbr;
}

MyPolygon *MyPolygon::gen_box(box &pix){
	return gen_box(pix.low[0],pix.low[1],pix.high[0],pix.high[1]);
}

vector<Point> MyPolygon::generate_test_points(int num){
	getMBB();
	vector<Point> points;
	for(int i=0;i<num;i++){
		Point p;
		p.x = get_rand_double()*(mbr->high[0]-mbr->low[0])+mbr->low[0];
		p.y = get_rand_double()*(mbr->high[1]-mbr->low[1])+mbr->low[1];
		points.push_back(p);
	}
	return points;
}


vector<MyPolygon *> MyPolygon::generate_test_polygons(int num){
	getMBB();
	vector<MyPolygon *> polys;
	for(int i=0;i<num;i++){
		double start_x = get_rand_double()*(mbr->high[0]-mbr->low[0])+mbr->low[0];
		double start_y = get_rand_double()*(mbr->high[1]-mbr->low[1])+mbr->low[1];

		double end_x = get_rand_double()*(mbr->high[0]-mbr->low[0])/100+start_x;
		double end_y = get_rand_double()*(mbr->high[1]-mbr->low[1])/100+start_y;
		MyPolygon *m = gen_box(start_x, start_y, end_x, end_y);
		polys.push_back(m);
	}

	return polys;
}

bool MyPolygon::convert_to_geos(geos::io::WKTReader *wkt_reader){
	try{
		//log("processing %d", ctx->source_polygons[i]->getid());
		geos_geom = wkt_reader->read(to_string(false, true));
		geos_geom->normalize();
		//log("%.10f",geos_geom->getArea());
	}catch(...){
		log("failed to parse polygon %ld", getid());
		print();
		geos_geom = NULL;
		return false;
	}
	delete boundary;
	boundary = NULL;

	return true;
}


void MyPolygon::print_without_head(bool print_hole, bool complete_ring){
	assert(boundary);
	cout<<"(";

	boundary->print(complete_ring);
	if(print_hole){
		for(VertexSequence *vs:this->holes){
			cout<<", ";
			vs->print(complete_ring);
		}
	}
	cout<<")";
}

void MyPolygon::print_triangles(){
	printf("MULTIPOLYGON(");
	bool first = true;
	for(int i=0;i<triangle_num;i++){
		if(!first){
			printf(",");
		}else{
			first = false;
		}
		Point *ps = triangles+3*i;
		printf("((%f %f, %f %f, %f %f))",ps[0].x,ps[0].y,ps[1].x,ps[1].y,ps[2].x,ps[2].y);
	}
	printf(")\n");
}

void MyPolygon::print(bool print_id, bool print_hole){
	if(print_id){
		cout<<"id:\t"<<this->id<<endl;
	}
	cout<<"POLYGON";
	print_without_head(print_hole);
	cout<<endl;
}
void MyPolygon::print_without_return(bool print_hole, bool complete_ring){
	cout<<"POLYGON";
	print_without_head(print_hole, complete_ring);
}

string MyPolygon::to_string(bool clockwise, bool complete_ring){
	std::stringstream ss;
	char double_str[200];
	ss<<"POLYGON";
	ss<<"((";
	if(clockwise){
		for(int i=0;i<boundary->num_vertices;i++){
			if(i!=0){
				ss<<",";
			}
			sprintf(double_str,"%f %f",boundary->p[i].x,boundary->p[i].y);
			ss<<double_str;
		}
		if(complete_ring){
			if(boundary->p[0].x!=boundary->p[boundary->num_vertices-1].x||
			   boundary->p[0].y!=boundary->p[boundary->num_vertices-1].y){
				ss<<",";
				sprintf(double_str,"%f %f",boundary->p[0].x,boundary->p[0].y);
				ss<<double_str;
			}
		}
	}else{
		for(int i=boundary->num_vertices-1;i>=0;i--){
			if(i!=boundary->num_vertices-1){
				ss<<",";
			}
			sprintf(double_str,"%f %f",boundary->p[i].x,boundary->p[i].y);
			ss<<double_str;
		}
		if(complete_ring){
			if(boundary->p[0].x!=boundary->p[boundary->num_vertices-1].x||
			   boundary->p[0].y!=boundary->p[boundary->num_vertices-1].y){
				ss<<",";
				sprintf(double_str,"%f %f",boundary->p[boundary->num_vertices-1].x,boundary->p[boundary->num_vertices-1].y);
				ss<<double_str;
			}
		}
	}

	ss<<"))";
	return ss.str();
}


/**
 *
 * 	if(qt.use_qtree){
		partition_qtree(qt.vpr);
		log("leaf count %d",qtree->leaf_count());
		log("size in bytes %d",qtree->size());
		std::stack<QTNode *> ws;
		ws.push(qtree);
		while(!ws.empty()){
			QTNode *cur = ws.top();
			ws.pop();
			if(cur->isleaf()){
				MyPolygon *m = gen_box(cur->mbr);
				if(cur->interior){
					inpolys->insert_polygon(m);
				}else if(cur->exterior){
					outpolys->insert_polygon(m);
				}else{
					borderpolys->insert_polygon(m);
				}
			}else{
				cur->push(ws);
			}
		}
	}
 *
 *
 */


void MyPolygon::print_partition(){
	MyMultiPolygon *inpolys = new MyMultiPolygon();
	MyMultiPolygon *borderpolys = new MyMultiPolygon();
	MyMultiPolygon *outpolys = new MyMultiPolygon();

	cout<<"border:"<<endl;
	borderpolys->print();
	cout<<"in:"<<endl;
	inpolys->print();
	cout<<"out:"<<endl;
	outpolys->print();

	delete borderpolys;
	delete inpolys;
	delete outpolys;
}

void MyPolygon::rasterization(int vpr){
	assert(vpr>0);
	if(raster){
		return;
	}
	pthread_mutex_lock(&ideal_partition_lock);
	if(raster==NULL){
		raster = new MyRaster(boundary,vpr);
		raster->rasterization();
	}
	pthread_mutex_unlock(&ideal_partition_lock);
}


QTNode *MyPolygon::partition_qtree(const int vpr){

	pthread_mutex_lock(&qtree_partition_lock);
	if(qtree){
		pthread_mutex_unlock(&qtree_partition_lock);
		return qtree;
	}

	int num_boxes = get_num_vertices()/vpr;
	int level = 1;
	while(pow(4,level)<num_boxes){
		level++;
	}
	int dimx = pow(2,level);
	int dimy = dimx;

	MyRaster *ras = new MyRaster(boundary,dimx,dimy);
	ras->rasterization();
	int box_count = 4;
	int cur_level = 1;

	qtree = new QTNode(*(this->getMBB()));
	std::stack<QTNode *> ws;
	qtree->split();
	qtree->push(ws);

	vector<QTNode *> level_nodes;

	//breadth first traverse
	query_context qc;
	while(box_count<num_boxes||!ws.empty()){
		if(cur_level>level){
			dimx *= 2;
			dimy *= 2;
			level = cur_level;
			delete ras;
			ras = new MyRaster(boundary,dimx,dimy);
			ras->rasterization();
		}

		while(!ws.empty()){
			QTNode *cur = ws.top();
			ws.pop();
			bool contained = false;
			if(!ras->contain(&cur->mbr, contained)){
				level_nodes.push_back(cur);
			}else if(contained){
				cur->interior = true;
			}else{
				cur->exterior = true;
			}
		}


		for(QTNode *n:level_nodes){
			if(box_count<num_boxes){
				n->split();
				n->push(ws);
				box_count += 3;
			}else{
				break;
			}
		}
		cur_level++;
		level_nodes.clear();
	}
	//assert(box_count>=num_boxes);

	delete ras;
	pthread_mutex_unlock(&qtree_partition_lock);

	return qtree;
}

size_t MyPolygon::raster_size(){
	size_t size = 0;
	const int nump = raster->get_num_pixels();

	const int nump_border = raster->get_num_pixels(BORDER);
	int bits_b = 0;
	int tmp = nump_border+2;
	while(tmp>0){
		bits_b++;
		tmp /= 2;
	}

	int numv = this->get_num_vertices();
	tmp = numv;
	int bits_v = 0;
	while(tmp>0){
		bits_v++;
		tmp /= 2;
	}

	bits_v = bits_v*2;

	int numc = raster->get_num_crosses();
	return (bits_b*nump + bits_v*nump_border + numc*64+7)/8;
}



