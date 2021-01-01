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
#include "../util/util.h"
#include "../index/RTree.h"


#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace std;

const static char *multipolygon_char = "MULTIPOLYGON";
const static char *polygon_char = "POLYGON";
const static char *point_char = "POINT";


enum PartitionStatus{
	OUT = 0,
	BORDER = 1,
	IN = 2
};

const static char *direction_str = "lrtb";

enum QueryType{
    contain = 0,
    distance = 1,
    within = 2
};

enum Direction{
	LEFT = 0,
	RIGHT = 1,
	TOP = 2,
	BOTTOM
};

class MyPolygon;

enum cross_type{
	ENTER = 0,
	LEAVE = 1
};

class cross_info{
public:
	cross_type type;
	double vertex;
};

class Point{
public:
	double x;
	double y;
	Point(){
		x = 0;
		y = 0;
	}
	Point(double xx, double yy){
		x = xx;
		y = yy;
	}
	void print(){
		printf("POINT (%f %f)\n",x,y);
		//printf("x: %f y: %f\n",x,y);
	}
	void print_without_return(){
		printf("POINT (%f %f)",x,y);
		//printf("x: %f y: %f\n",x,y);
	}
	string to_string(){
		std::stringstream ss;
		char double_str[200];
		ss<<"POINT("<<x<<" "<<y<<")";
		return ss.str();
	}
	static Point *read_one_point(string &str);
};

class Pixel{
public:
	int id[2];
	double low[2] = {10000000,10000000};
	double high[2] = {-10000000,-10000000};
	PartitionStatus status = OUT;

	int vstart = -1;
	int vend = -1;
	Pixel(){

	}

	void update(Point &p){
		if(low[0]>p.x){
			low[0] = p.x;
		}
		if(high[0]<p.x){
			high[0] = p.x;
		}

		if(low[1]>p.y){
			low[1] = p.y;
		}
		if(high[1]<p.y){
			high[1] = p.y;
		}
	}

	bool intersect(Pixel *target){
		return !(target->low[0]>high[0]||
				 target->high[0]<low[0]||
				 target->low[1]>high[1]||
				 target->high[1]<low[1]);
	}
	bool contain(Pixel *target){
		return target->low[0]>=low[0]&&
			   target->high[0]<=high[0]&&
			   target->low[1]>=low[1]&&
			   target->high[1]<=high[1];
	}
	bool contain(Point &p){
		return p.x>=low[0]&&
			   p.x<=high[0]&&
			   p.y>=low[1]&&
			   p.y<=high[1];
	}

	double max_distance(Point p){
		double md = 0;
		double dist = (p.x-low[0])*(p.x-low[0])+(p.y-low[1])*(p.y-low[1]);
		if(dist>md){
			md = dist;
		}
		dist = (p.x-low[0])*(p.x-low[0])+(p.y-high[1])*(p.y-high[1]);
		if(dist>md){
			md = dist;
		}
		dist = (p.x-high[0])*(p.x-high[0])+(p.y-low[1])*(p.y-low[1]);
		if(dist>md){
			md = dist;
		}
		dist = (p.x-high[0])*(p.x-high[0])+(p.y-high[1])*(p.y-high[1]);
		if(dist>md){
			md = dist;
		}
		return sqrt(md);
	}

	double area(){
		return (high[0]-low[0])*(high[1]-low[1]);
	}

	double distance(Point &p){
		if(this->contain(p)){
			return 0;
		}
		double dx = max(abs(p.x-(low[0]+high[0])/2) - (high[0]-low[0])/2, 0.0);
		double dy = max(abs(p.y-(low[1]+high[1])/2) - (high[1]-low[1])/2, 0.0);
		return sqrt(dx * dx + dy * dy);
	}

    double distance_geography(Point &p){
    	if(this->contain(p)){
    		return 0.0;
    	}
        double dx = max(abs(p.x-(low[0]+high[0])/2) - (high[0]-low[0])/2, 0.0);
        double dy = max(abs(p.y-(low[1]+high[1])/2) - (high[1]-low[1])/2, 0.0);
        dx = dx/degree_per_kilometer_latitude;
        dy = dy/degree_per_kilometer_longitude(p.y);
        return sqrt(dx * dx + dy * dy);
    }

	vector<cross_info> crosses[4];
	void enter(double val, Direction d, int vnum);
	void leave(double val, Direction d, int vnum);

	MyPolygon *to_polygon();

	void print(){
		switch(status){
		case IN:
			printf("I-");
			break;
		case OUT:
			printf("O-");
			break;
		case BORDER:
			printf("B-");
			break;
		}
		printf("%3d %3d-",id[0],id[1]);
		printf("low_x: %f low_y %f high_x %f high_y %f\n",low[0],low[1],high[0],high[1]);
	}

};




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
	VertexSequence *clone();
	void print();
	bool clockwise();
	void reverse();
	double area();
	bool contain(Point p);
	//fix the ring
	void fix();
};



enum QT_Direction{
	bottom_left = 0,
	bottom_right = 1,
	top_left = 2,
	top_right = 3
};

class QTNode{

public:

	bool isleaf = true;
	QTNode *children[4];
	Pixel mbr;
	bool interior = false;
	bool exterior = false;
	int level = 0;

	QTNode(double low_x, double low_y, double high_x, double high_y);
	QTNode(Pixel m);
	~QTNode();

	void split();
	void push(std::stack<QTNode *> &ws);
	int size();
	int leaf_count();

	bool within(Point &p, double within_dist);
	void retrieve_within(Point &p, vector<QTNode *> &candidates, double within_dist);
	QTNode *retrieve(Point &p);
	bool determine_contain(Point &p);
	bool determine_contain(Pixel &p, bool &has_in, bool &has_ex);
	bool determine_contain(Pixel &p);

};

class query_context{
public:
	//configuration
	int thread_id = 0;
	int num_threads = 0;
	int vpr = 10;
	int vpr_end = 10;
	bool use_grid = false;
	bool use_qtree = false;
	bool use_convex_hull = false;
	bool use_mer = false;
	bool use_triangulate = false;
	int mer_sample_round = 20;
	bool perform_refine = true;
	bool gpu = false;
	bool collect_latency = false;
	float sample_rate = 1.0;
	int small_threshold = 500;
	int big_threshold = 1000000;
	bool sort_polygons = false;
	int report_gap = 1000000;
	int distance_buffer_size = 10;
	QueryType query_type = QueryType::contain;
	string source_path;
	string target_path;

	string valid_path;

	//query target, for temporary use
	void *target = NULL;
	void *target2 = NULL;
	query_context *global_ctx = NULL;

	//shared staff
	size_t index = 0;
	size_t index_end = 0;
	pthread_mutex_t lock;


	//result
	double distance = 0;
	bool filter_checked_only = false;

	//query statistic
	size_t found = 0;
	size_t query_count = 0;
	size_t checked_count = 0;
	size_t refine_count = 0;

	size_t pixel_checked = 0;
	size_t border_checked = 0;
	size_t edges_checked = 0;
	double pixel_check_time = 0;
	double edges_check_time = 0;
	double check_time = 0;


	// temporary storage for query processing
	vector<MyPolygon *> source_polygons;
	vector<MyPolygon *> target_polygons;
	double *points = NULL;
	size_t target_num = 0;


	vector<pair<MyPolygon *, MyPolygon *>> candidates;
	map<int, int> vertex_number;
	map<int, double> latency;


	// functions
	query_context();
	~query_context();
	query_context(query_context &t);
	query_context operator + (query_context const &obj);
	query_context& operator=(query_context const &obj);

	void report_latency(int num_v, double latency);
	void load_points();
	void report_progress();
	void merge_global();
	bool next_batch(int batch_num=1);
	void reset_stats(){
		//query statistic
		found = 0;
		query_count = 0;
		checked_count = 0;
		refine_count = 0;
		pixel_checked = 0;
		border_checked = 0;
		edges_checked = 0;
		pixel_check_time = 0;
		edges_check_time = 0;
		check_time = 0;
		index = 0;
	}

	bool is_within_query(){
		return query_type == QueryType::within;
	}

	bool is_contain_query(){
		return query_type == QueryType::contain;
	}
	void print_stats();
};


query_context get_parameters(int argc, char **argv);
int *polygon_triangulate ( int n, double x[], double y[] );

class MyPolygon{
	Pixel *mbr = NULL;
	Pixel *mer = NULL;
	vector<vector<Pixel>> partitions;
	QTNode *qtree = NULL;
	RTree<int *, double, 2, double> *rtree = NULL;
	int *triangles = NULL;

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

	int *triangulate();
	void build_rtree();



	static VertexSequence *read_vertices(const char *wkt, size_t &offset, bool clockwise=true);
	static MyPolygon *read_polygon(const char *wkt, size_t &offset);
	static MyPolygon *gen_box(double minx,double miny,double maxx,double maxy);
	static MyPolygon *read_one_polygon();

	static vector<MyPolygon *> load_binary_file(const char *path, query_context ctx, bool sample=false);
	static MyPolygon * load_binary_file_single(const char *path, query_context ctx, int idx);
	static MyPolygon * read_polygon_binary_file(ifstream &is);

	bool contain(Point p);// brute-forcely check containment
	bool contain(Point p, query_context *ctx);
	bool contain_rtree(Point p, query_context *ctx);
	bool intersect(MyPolygon *target, query_context *ctx);
	bool intersect_segment(Pixel *target);
	bool contain(MyPolygon *target, query_context *ctx);
	bool contain(Pixel *target, query_context *ctx);
	double distance(Point &p, query_context *ctx);
	double distance(Point &p);// brute-forcely calculate distance

	bool contain_try_partition(Pixel *target, query_context *ctx);

	void print_without_head(bool print_hole = false);
	void print(bool print_hole=false);
	void print_without_return(bool print_hole=false);
	string to_string(bool clockwise = true);
	Pixel *getMBB();
	Pixel *getMER(){
		query_context ctx;
		return getMER(&ctx);
	}
	Pixel *getMER(query_context *ctx);
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
		for(VertexSequence *vs:internal_polygons){
			double a = vs->area();
			//assert(a<=0);
			sum += a;
		}
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


#endif /* SRC_MYPOLYGON_H_ */
