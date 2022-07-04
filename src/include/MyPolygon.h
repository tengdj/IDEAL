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

#include "../geometry/triangulate/poly2tri.h"
#include "util.h"
#include "../index/RTree.h"
#include "../index/QTree.h"
#include "Pixel.h"
#include "Point.h"
#include "query_context.h"

using namespace std;
using namespace p2t;

const static char *multipolygon_char = "MULTIPOLYGON";
const static char *polygon_char = "POLYGON";

class VertexSequence{
public:
	int num_vertices = 0;
	Point *p = NULL;
public:

	VertexSequence(){};
	VertexSequence(int nv);
	VertexSequence(int nv, Point *pp);
	size_t encode(char *dest);
	size_t decode(char *source);
	size_t get_data_size();

	~VertexSequence();
	vector<Vertex *> pack_to_polyline();
	VertexSequence *clone();
	VertexSequence *convexHull();
	box *getMBR();
	void print(bool complete_ring=false);
	bool clockwise();
	void reverse();
	double area();
	bool contain(Point &p);
	//fix the collinear problem
	void fix();
};


class MyRaster{
	box *mbr = NULL;
	VertexSequence *vs = NULL;
	vector<vector<Pixel *>> pixels;
	double step_x = 0.0;
	double step_y = 0.0;
	int dimx = 0;
	int dimy = 0;
	void init_pixels();
	void evaluate_edges();
	void scanline_reandering();

public:

	MyRaster(VertexSequence *vs, int epp);
	MyRaster(VertexSequence *vs, int dimx, int dimy);
	void rasterization();
	~MyRaster();

	bool contain(box *,bool &contained);
	vector<Pixel *> get_intersect_pixels(box *pix);
	vector<Pixel *> get_closest_pixels(box *target);
	Pixel *get_pixel(Point &p);
	Pixel *get_closest_pixel(Point &p);
	Pixel *get_closest_pixel(box *target);
	vector<Pixel *> expand_radius(int lowx, int highx, int lowy, int highy, int step);
	vector<Pixel *> expand_radius(Pixel *center, int step);

	int get_offset_x(double x);
	int get_offset_y(double y);

	/* statistics collection*/
	int count_intersection_nodes(Point &p);
	int get_num_border_edge();
	size_t get_num_pixels();
	size_t get_num_pixels(PartitionStatus status);
	size_t get_num_gridlines();
	size_t get_num_crosses();
	void print();

	vector<Pixel *> get_pixels(PartitionStatus status);
	box *extractMER(Pixel *starter);

	vector<Pixel *> retrieve_pixels(box *);

	/*
	 * the gets functions
	 *
	 * */
	double get_step_x(){
		return step_x;
	}
	double get_step_y(){
		return step_y;
	}
	double get_step(bool geography){
		if(geography){
			return min(step_x/degree_per_kilometer_longitude(mbr->low[1]), step_y/degree_per_kilometer_latitude);
		}else{
			return min(step_x, step_y);
		}
	}
	int get_dimx(){
		return dimx;
	}
	int get_dimy(){
		return dimy;
	}

	Pixel *get(int dx, int dy){
		assert(dx>=0&&dx<=dimx);
		assert(dy>=0&&dy<=dimy);
		return pixels[dx][dy];
	}

};

// the structured metadata of a polygon
typedef struct PolygonMeta_{
	uint size; // size of the polygon in bytes
	uint num_vertices; // number of vertices in the boundary (hole excluded)
	size_t offset; // the offset in the file
	box mbr; // the bounding boxes
} PolygonMeta;

class MyPolygon{
	size_t id = 0;

	box *mbr = NULL;
	box *mer = NULL;
	MyRaster *raster = NULL;

	QTNode *qtree = NULL;
	// for triangulation
	Point *triangles = NULL;
	size_t triangle_num = 0;


	RTNode *rtree = NULL;

	pthread_mutex_t ideal_partition_lock;
	pthread_mutex_t qtree_partition_lock;


public:
	// the Hilbert curve value of the MBR centroid
	size_t hc_id = 0;
	VertexSequence *boundary = NULL;
	VertexSequence *convex_hull = NULL;
	vector<VertexSequence *> holes;
	MyPolygon(){
	    pthread_mutex_init(&ideal_partition_lock, NULL);
	    pthread_mutex_init(&qtree_partition_lock, NULL);
	}
	~MyPolygon();
	MyPolygon *clone();
	MyRaster *get_rastor(){
		return raster;
	}
	size_t get_num_pixels(){
		if(raster){
			return raster->get_num_pixels();
		}
		return 0;
	}
	size_t get_num_pixels(PartitionStatus status){
		if(raster){
			return raster->get_num_pixels(status);
		}
		return 0;
	}
	double get_pixel_portion(PartitionStatus status){
		if(raster){
			return 1.0*get_num_pixels(status)/get_num_pixels();
		}
		return 0.0;
	}

	void triangulate();
	size_t get_triangle_size(){
		size_t sz = 0;
		size_t num_bits = 0;
		int numv = this->get_num_vertices();
		while(numv>0){
			num_bits++;
			numv/=2;
		}
		sz += triangle_num*3*(num_bits+7)/8;
		return sz;
	}
	void build_rtree();
	RTNode *get_rtree(){
		return rtree;
	}
	size_t get_rtree_size(){
		if(rtree){
			int count = rtree->node_count();
			return count*4*8;
		}else{
			return 0;
		}
	}


	static VertexSequence *read_vertices(const char *wkt, size_t &offset, bool clockwise=true);
	static MyPolygon *read_polygon(const char *wkt, size_t &offset);
	static bool validate_polygon(const char *wkt, size_t &offset, size_t &len);
	static bool validate_vertices(const char *wkt, size_t &offset, size_t &len);

	static MyPolygon *gen_box(double minx,double miny,double maxx,double maxy);
	static MyPolygon *gen_box(box &pix);
	static MyPolygon *read_one_polygon();

	bool contain(Point &p);// brute-forcely check containment
	bool contain(Point &p, query_context *ctx, bool profile = true);
	bool intersect(MyPolygon *target, query_context *ctx);
	bool intersect_segment(box *target);
	bool contain(MyPolygon *target, query_context *ctx);
	double distance(Point &p, query_context *ctx, bool profile = true);
	double distance(MyPolygon *target, query_context *ctx);
	double distance(MyPolygon *target, Pixel *pix, query_context *ctx, bool profile = true);

	double distance_rtree(Point &p, query_context *ctx);
	double distance_rtree(Point &start, Point &end, query_context *ctx);
	bool intersect_rtree(Point &start, Point &end, query_context *ctx);


	void print_without_head(bool print_hole = false, bool complete_ring = false);
	void print(bool print_id=true, bool print_hole=false);
	void print_triangles();
	void print_without_return(bool print_hole=false, bool complete_ring=false);
	string to_string(bool clockwise = false, bool complete_ring=false);
	box *getMBB();
	box *getMER(query_context *ctx=NULL);

	VertexSequence *get_convex_hull();

	size_t partition_size();

	void rasterization(int vertex_per_raster);

	QTNode *partition_qtree(const int vpr);
	QTNode *get_qtree(){
		return qtree;
	}

	inline int get_num_vertices(){
		if(!boundary){
			return 0;
		}
		return boundary->num_vertices;
	}
	inline double getx(int index){
		assert(boundary&&index<boundary->num_vertices);
		return boundary->p[index].x;
	}
	inline double gety(int index){
		assert(boundary&&index<boundary->num_vertices);
		return boundary->p[index].y;
	}
	inline Point *get_point(int index){
		assert(boundary&&index<boundary->num_vertices);
		return &boundary->p[index];
	}

	inline size_t getid(){
		return id;
	}
	inline void setid(size_t iid){
		id = iid;
	}

	PolygonMeta get_meta();
	size_t get_data_size();
	size_t encode(char *target);
	size_t decode(char *source);
	static char *encode_raster(vector<vector<Pixel>> raster);
	static vector<vector<Pixel>> decode_raster(char *);

	vector<Point> generate_test_points(int num);
	vector<MyPolygon *> generate_test_polygons(int num);

	double area(){
		double area_buffer = 0;
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

class MyMultiPolygon{
	vector<MyPolygon *> polygons;
public:
	static MyMultiPolygon *read_multipolygon();
	static MyPolygon *read_one_polygon();
	static bool validate_wkt(string &wkt_str);

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

class MyMultiPoint{
	vector<Point *> points;
public:
	void insert(vector<Point *> &origin_points);
	void insert(Point *p);
	void print(size_t max_num = INT_MAX);
};


//utility functions
void process_rasterization(query_context *ctx);
void process_convex_hull(query_context *ctx);
void process_mer(query_context *ctx);
void process_internal_rtree(query_context *gctx);
void preprocess(query_context *gctx);

void print_boxes(vector<box *> boxes);


// storage related functions
size_t load_points_from_path(const char *path, Point **points);
size_t load_mbr_from_file(const char *path, box **);
size_t load_polygonmeta_from_file(const char *path, PolygonMeta **pmeta);

void dump_to_file(const char *path, char *data, size_t size);
void dump_polygons_to_file(vector<MyPolygon *> polygons, const char *path);
vector<MyPolygon *> load_binary_file(const char *path, query_context &ctx);
size_t number_of_objects(const char *path);
box universe_space(const char *path);

#endif /* SRC_MYPOLYGON_H_ */
