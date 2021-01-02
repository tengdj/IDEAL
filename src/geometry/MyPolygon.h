/*
 * MyPolygon.h
 *
 *  Created on: May 9, 2020
 *      Author: teng
 */

#ifndef SRC_MYPOLYGON_H_
#define SRC_MYPOLYGON_H_

#include <vector>
#include <string>
#include <stdint.h>
#include <math.h>
#include <stack>
#include <map>
#include <bits/stdc++.h>
#include <boost/program_options.hpp>

#include "../triangulate/poly2tri.h"
#include "../util/util.h"
#include "../index/RTree.h"
#include "../index/QTree.h"
#include "Pixel.h"
#include "Point.h"
#include "query_context.h"

namespace po = boost::program_options;
using namespace std;
using namespace p2t;

const static char *multipolygon_char = "MULTIPOLYGON";
const static char *polygon_char = "POLYGON";

class VertexSequence{
public:
	int num_vertices = 0;
	double *x = NULL;
	double *y = NULL;
	VertexSequence(){};
	VertexSequence(int nv);
	VertexSequence(int nv, double *xx, double *yy);
	size_t encode(char *dest);
	size_t decode(char *source);
	size_t get_data_size();

	~VertexSequence();
	vector<Point *> pack_to_polyline();
	VertexSequence *clone();
	void print();
	bool clockwise();
	void reverse();
	double area();
	bool contain(Point &p);
	//fix the ring
	void fix();
};




class MyPolygon{
	Pixel *mbr = NULL;
	Pixel *mer = NULL;
	vector<vector<Pixel>> partitions;
	QTNode *qtree = NULL;
	RTree<Triangle *, double, 2, double> *rtree = NULL;
	Pixel *rtree_pixel = NULL;
	CDT *cdt = NULL;;

	double step_x = 0;
	double step_y = 0;
	int id = 0;
	double area_buffer = -1;
	pthread_mutex_t ideal_partition_lock;
	pthread_mutex_t qtree_partition_lock;


public:
	size_t offset = 0;
	VertexSequence *boundary = NULL;
	VertexSequence *convex_hull = NULL;
	bool valid_for_triangulate = false;
	vector<VertexSequence *> internal_polygons;

	MyPolygon(){
	    pthread_mutex_init(&ideal_partition_lock, NULL);
	    pthread_mutex_init(&qtree_partition_lock, NULL);
	}
	~MyPolygon();
	MyPolygon *clone();

	void triangulate();
	size_t get_triangle_size(){
		size_t sz = 0;
		size_t num_bits = 0;
		int numv = this->get_num_vertices();
		while(numv>0){
			num_bits++;
			numv/=2;
		}
		if(cdt){
			sz += cdt->GetTriangles().size()*3*(num_bits+7)/8;
		}
		return sz;
	}
	RTree<Triangle *, double, 2, double> * build_rtree();
	RTree<Triangle *, double, 2, double> * get_rtree(){
		return rtree;
	}
	Pixel *get_rtree_pixel(){
		return this->rtree_pixel;
	}
	size_t get_rtree_size(){
		if(rtree){
			return rtree->Count()*4*8;
		}else{
			return 0;
		}
	}


	static VertexSequence *read_vertices(const char *wkt, size_t &offset, bool clockwise=true);
	static MyPolygon *read_polygon(const char *wkt, size_t &offset);
	static MyPolygon *gen_box(double minx,double miny,double maxx,double maxy);
	static MyPolygon *gen_box(Pixel &pix);
	static MyPolygon *read_one_polygon();

	static vector<MyPolygon *> load_binary_file(const char *path, query_context ctx, bool sample=false);
	static MyPolygon * load_binary_file_single(const char *path, query_context ctx, int idx);
	static MyPolygon * read_polygon_binary_file(ifstream &is);

	bool contain(Point &p);// brute-forcely check containment
	bool contain(Point &p, query_context *ctx);
	bool contain_rtree(Point &p, query_context *ctx);
	bool intersect(MyPolygon *target, query_context *ctx);
	bool intersect_segment(Pixel *target);
	bool contain(MyPolygon *target, query_context *ctx);
	bool contain(Pixel *target, query_context *ctx);
	double distance(Point &p, query_context *ctx);
	double distance(Point &p);// brute-forcely calculate distance
	double distance_rtree(Point &p, query_context *ctx);

	bool contain_try_partition(Pixel *target, query_context *ctx);

	void print_without_head(bool print_hole = false);
	void print(bool print_hole=false);
	void print_triangles();
	void print_without_return(bool print_hole=false);
	string to_string(bool clockwise = true);
	Pixel *getMBB();
	Pixel *getMER(query_context *ctx=NULL);
	Pixel *generateMER(int cx, int cy);

	void init_partition(const int dimx, const int dimy);
	void evaluate_edges(const int dimx, const int dimy);


	VertexSequence *get_convex_hull();

	size_t partition_size();
	size_t get_cross_num();
	size_t get_grid_num(){
		if(partitions.size()==0){
			return 0;
		}
		return partitions.size()+partitions[0].size();
	}
	vector<vector<Pixel>> partition(int vertex_per_raster);
	vector<vector<Pixel>> partition(int xdim, int ydim);
	vector<vector<Pixel>> partition_with_query(int vertex_per_raster);
	vector<Pixel *> get_pixels(PartitionStatus status);

	Pixel *get_closest_pixel(Point p);
	QTNode *partition_qtree(const int vpr);
	QTNode *get_qtree(){
		return qtree;
	}

	void reset_grid_partition(){
		for(vector<Pixel> &rows:partitions){
			rows.clear();
		}
		partitions.clear();
	}

	bool is_grid_partitioned(){
		return partitions.size()>0;
	}
	bool is_qtree_partitioned(){
		return qtree!=NULL;
	}
	inline int get_num_vertices(){
		if(!boundary){
			return 0;
		}
		return boundary->num_vertices;
	}
	inline double getx(int index){
		assert(boundary&&index<boundary->num_vertices);
		return boundary->x[index];
	}
	inline double gety(int index){
		assert(boundary&&index<boundary->num_vertices);
		return boundary->y[index];
	}
	inline int get_pixel_x(double xval){
		assert(mbr);
		int x = double_to_int((xval-mbr->low[0])/step_x);
		return x;
	}
	inline int get_pixel_y(double yval){
		assert(mbr);
		int y = double_to_int((yval-mbr->low[1])/step_y);
		return y;
	}
	inline int getid(){
		return id;
	}
	inline void setid(int iid){
		id = iid;
	}

	size_t get_data_size();
	size_t encode_to(char *target);
	size_t decode_from(char *source);

	int get_num_partitions(){
		if(partitions.size()==0){
			return 0;
		}
		return partitions.size()*partitions[0].size();
	}
	int get_num_partitions(PartitionStatus status){
		if(partitions.size()==0){
			return 0;
		}
		int num = 0;
		for(vector<Pixel> &rows:partitions){
			for(Pixel &p:rows){
				if(p.status==status){
					num++;
				}
			}
		}
		return num;
	}

	static char *encode_partition(vector<vector<Pixel>> partitions);
	static vector<vector<Pixel>> decode_partition(char *);
	vector<Point> generate_test_points(int num);
	vector<MyPolygon *> generate_test_polygons(int num);

	double area(){
		if(area_buffer>=0){
			return area_buffer;
		}
		assert(boundary);
		double sum = boundary->area();
		//assert(sum>=0);
//		for(VertexSequence *vs:internal_polygons){
//			double a = vs->area();
//			//assert(a<=0);
//			sum += a;
//		}
		//assert(sum>=0);
		area_buffer = sum;
		return area_buffer/2;
	}

	void print_partition(query_context qc);
};

class queue_element{
public:
	int id = 0;
	vector<MyPolygon *> polys;
	queue_element(int i){
		id = i;
	}
	~queue_element(){
		for(MyPolygon *p:polys){
			delete p;
		}
		polys.clear();
	}
};

class MyMultiPolygon{
	vector<MyPolygon *> polygons;
public:
	static MyPolygon * read_one_polygon();
	MyMultiPolygon(const char *wkt);
	MyMultiPolygon(){};
	~MyMultiPolygon();
	int num_polygons(){
		return polygons.size();
	}
	void print();
	void print_separate();
	vector<MyPolygon *> get_polygons(){
		return polygons;
	}
	MyPolygon *get_polygon(int index){
		return index>=polygons.size()?NULL:polygons[index];
	}
	size_t insert_polygon(MyPolygon *p){
		polygons.push_back(p);
		return polygons.size();
	}
	size_t insert_polygon(vector<MyPolygon *> ps){
		for(MyPolygon *p:ps){
			polygons.push_back(p);
		}
		return polygons.size();
	}
	size_t get_data_size(){
		size_t ds = 0;
		for(MyPolygon *p:polygons){
			ds += p->get_data_size();
		}
		return ds;
	}
};


//utility functions

void process_partition(query_context *ctx);
void process_convex_hull(query_context *ctx);
void process_mer(query_context *ctx);
void process_triangulate(query_context *gctx);
void process_internal_rtree(query_context *gctx);
void preprocess(query_context *gctx);
#endif /* SRC_MYPOLYGON_H_ */
