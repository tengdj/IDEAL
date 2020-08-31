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
	Direction direction;
	double vertex;
	static bool overwrites(cross_info &enter1, cross_info &leave1, cross_info &enter2, cross_info &leave2);
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
	static Point *read_one_point(string &str);
};

class Pixel{
public:
	int id[2];
	double low[2];
	double high[2];
	PartitionStatus status = OUT;
	PartitionStatus border[4];
	int vstart = -1;
	int vend = -1;
	Pixel(){
		border[0] = OUT;
		border[1] = OUT;
		border[2] = OUT;
		border[3] = OUT;
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

	double distance(Point &p){
		double dx = max(abs(p.x-(low[0]+high[0])/2) - (high[0]-low[0])/2, 0.0);
		double dy = max(abs(p.y-(low[1]+high[1])/2) - (high[1]-low[1])/2, 0.0);
		return sqrt(dx * dx + dy * dy);
	}

	vector<cross_info> crosses;
	void enter(double val, Direction d, int vnum);
	void leave(double val, Direction d, int vnum);
	void process_enter_leave();
	void setin(){
		status = IN;
		border[0] = IN;
		border[1] = IN;
		border[2] = IN;
		border[3] = IN;
	}

	MyPolygon *to_polygon();
	bool overwrites(cross_info &enter1, cross_info &leave1, cross_info &enter2, cross_info &leave2);
	void print(){
		printf("low_x: %f low_y %f high_x %f high_y %f\n",low[0],low[1],high[0],high[1]);
	}

};




class VertexSequence{
public:
	int num_vertices = 0;
	double *x = NULL;
	double *y = NULL;
	VertexSequence(){};
	VertexSequence(int nv){
		x = new double[nv];
		y = new double[nv];
		num_vertices = nv;
	}
	VertexSequence(int nv, double *xx, double *yy){
		num_vertices = nv;
		x = new double[nv];
		y = new double[nv];
		memcpy((char *)x,(char *)xx,num_vertices*sizeof(double));
		memcpy((char *)y,(char *)yy,num_vertices*sizeof(double));
	};


	size_t encode(char *dest){
		size_t encoded = 0;
		((long *)dest)[0] = num_vertices;
		encoded += sizeof(long);
		memcpy(dest+encoded,(char *)x,num_vertices*sizeof(double));
		encoded += num_vertices*sizeof(double);
		memcpy(dest+encoded,(char *)y,num_vertices*sizeof(double));
		encoded += num_vertices*sizeof(double);
		return encoded;
	}

	size_t decode(char *source){
		size_t decoded = 0;
		num_vertices = ((long *)source)[0];
		x = new double[num_vertices];
		y = new double[num_vertices];
		decoded += sizeof(long);
		memcpy((char *)x,source+decoded,num_vertices*sizeof(double));
		decoded += num_vertices*sizeof(double);
		memcpy((char *)y,source+decoded,num_vertices*sizeof(double));
		decoded += num_vertices*sizeof(double);
		return decoded;
	}

	size_t get_data_size(){
		return sizeof(long)+num_vertices*2*sizeof(double);
	}

	~VertexSequence(){
		if(x){
			delete []x;
		}
		if(y){
			delete []y;
		}
	}
	VertexSequence *clone(){
		VertexSequence *ret = new VertexSequence();
		ret->num_vertices = num_vertices;
		ret->x = new double[num_vertices];
		ret->y = new double[num_vertices];
		memcpy((void *)ret->x,(void *)x,sizeof(double)*num_vertices);
		memcpy((void *)ret->y,(void *)y,sizeof(double)*num_vertices);
		return ret;
	}
	bool clockwise();
	void reverse();
	double area();
};

class query_context{
public:
	int thread_id = 0;
	void *target = NULL;
	Point target_p;
	int vpr = 10;
	bool use_grid = false;
	bool use_qtree = false;
	bool query_vector = true;

	int found = 0;
	int checked = 0;
	double distance = 0;
	bool gpu = false;
	size_t query_count = 0;
	size_t inout_check = 0;
	size_t border_check = 1;
	size_t edges_checked = 0;
	bool rastor_only = false;
	float sample_rate = 1.0;

	int small_threshold = 500;
	int big_threshold = 1000000;

	size_t total_data_size = 0;
	size_t total_partition_size = 0;

	size_t partitions_count = 0;

	bool sort_polygons = false;

	vector<pair<MyPolygon *, MyPolygon *>> candidates;
	map<int, int> partition_vertex_number;
	map<int, double> partition_latency;

	query_context(){}
	~query_context(){candidates.clear();}
	query_context operator + (query_context const &obj) {
		query_context res = *this;
		res.query_count += obj.query_count;
		res.inout_check += obj.inout_check;
		res.border_check += obj.border_check;
		res.found += obj.found;
		res.edges_checked += obj.edges_checked;
		res.total_data_size += obj.total_data_size;
		res.total_partition_size += obj.total_partition_size;
		res.partitions_count += obj.partitions_count;
		for(auto &it :obj.partition_vertex_number){

			const double lt = obj.partition_latency.at(it.first);
			if(res.partition_vertex_number.find(it.first)!=res.partition_vertex_number.end()){
				res.partition_vertex_number[it.first] = res.partition_vertex_number[it.first]+it.second;
				res.partition_latency[it.first] = res.partition_latency[it.first]+lt;
			}else{
				res.partition_vertex_number[it.first] = it.second;
				res.partition_latency[it.first] = lt;
			}
		}
		return res;
	}
	void report_latency(int num_v, double latency){
		if(partition_vertex_number.find(num_v)==partition_vertex_number.end()){
			partition_vertex_number[num_v] = 1;
			partition_latency[num_v] = latency;
		}else{
			partition_vertex_number[num_v] = partition_vertex_number[num_v]+1;
			partition_latency[num_v] = partition_latency[num_v]+latency;
		}

	}
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
	Pixel mbb;
	bool interior = false;
	bool exterior = false;
	int level = 0;

	QTNode(double low_x, double low_y, double high_x, double high_y){
		mbb.low[0] = low_x;
		mbb.low[1] = low_y;
		mbb.high[0] = high_x;
		mbb.high[1] = high_y;
	}
	QTNode(Pixel m){
		mbb = m;
	}
	void split(){
		isleaf = false;
		double mid_x = (mbb.high[0]+mbb.low[0])/2;
		double mid_y = (mbb.high[1]+mbb.low[1])/2;
		children[bottom_left] = new QTNode(mbb.low[0],mbb.low[1],mid_x,mid_y);
		children[bottom_right] = new QTNode(mid_x,mbb.low[1],mbb.high[0],mid_y);
		children[top_left] = new QTNode(mbb.low[0],mid_y,mid_x,mbb.high[1]);
		children[top_right] = new QTNode(mid_x,mid_y,mbb.high[0],mbb.high[1]);
		for(int i=0;i<4;i++){
			children[i]->level = level+1;
		}
	}
	void push(std::stack<QTNode *> &ws){
		if(!isleaf){
			for(int i=0;i<4;i++){
				ws.push(children[i]);
			}
		}
	}
	~QTNode(){
		if(!isleaf){
			for(int i=0;i<4;i++){
				delete children[i];
			}
		}
	}
	int size(){
		int cs = 2+4*8;
		if(!isleaf){
			for(int i=0;i<4;i++){
				cs += children[i]->size();
			}
		}
		return cs;
	}
	int leaf_count(){
		if(isleaf){
			return 1;
		}else{
			int lc = 0;
			for(int i=0;i<4;i++){
				lc += children[i]->leaf_count();
			}
			return lc;
		}
	}
	QTNode *retrieve(Point &p){
		if(this->isleaf){
			return this;
		}else{
			int offset =  2*(p.y>(mbb.low[1]+mbb.high[1])/2)+(p.x>(mbb.low[0]+mbb.high[0])/2);
			return children[offset]->retrieve(p);
		}
	}

	bool determine_contain(Point &p){
		assert(mbb.contain(p));
		QTNode *n = retrieve(p);
		return n->interior||n->exterior;
	}

	bool determine_contain(Pixel &p, bool &has_in, bool &has_ex){
		if(!mbb.intersect(&p)){
			return true;
		}
		if(isleaf){
			//border box
			if(!interior&&!exterior){
				return false;
			}
			has_ex |= exterior;
			has_in |= interior;
			//cross box
			if(has_in&&has_ex){
				return false;
			}
		}else{
			for(int i=0;i<4;i++){
				if(!children[i]->determine_contain(p, has_in, has_ex)){
					return false;
				}
			}
		}
		return true;
	}

	bool determine_contain(Pixel &p){
		assert(this->level==0);
		bool has_in = false;
		bool has_out = false;
		return determine_contain(p, has_in, has_out);
	}


};

class MyPolygon{
	Pixel *mbb = NULL;
	vector<vector<Pixel>> partitions;
	QTNode *qtree = NULL;
	double step_x = 0;
	double step_y = 0;
	int id = 0;
	double area_buffer = -1;
	pthread_mutex_t partition_lock;

public:
	unsigned int offset = 0;
	VertexSequence *boundary = NULL;
	vector<VertexSequence *> internal_polygons;


	static VertexSequence *read_vertices(const char *wkt, size_t &offset, bool clockwise=true);
	static MyPolygon *read_polygon(const char *wkt, size_t &offset);
	static MyPolygon *gen_box(double minx,double miny,double maxx,double maxy);
	static MyPolygon *read_one_polygon();

	static vector<MyPolygon *> load_binary_file(const char *path, query_context ctx);
	static MyPolygon * read_polygon_binary_file(ifstream &is);

	bool contain(Point p, query_context *ctx);
	bool intersect(MyPolygon *target, query_context *ctx);
	bool intersect_segment(Pixel *target);
	bool contain(MyPolygon *target, query_context *ctx);
	bool contain(Pixel *target, query_context *ctx);
	double distance(Point &p, query_context *ctx);

	bool contain_try_partition(Pixel *target, query_context *ctx);

	void internal_partition();
	MyPolygon *clone();
	MyPolygon(){
	    pthread_mutex_init(&partition_lock, NULL);
	}
	~MyPolygon();
	void print_without_head(bool print_hole = false);
	void print(bool print_hole=false);
	void print_without_return(bool print_hole=false);
	string to_string(bool clockwise = true);
	Pixel *getMBB();

	void evaluate_border(const int dimx, const int dimy);
	void spread_pixels(const int dimx, const int dimy);

	vector<vector<Pixel>> partition(int vertex_per_raster);
	vector<vector<Pixel>> partition(int xdim, int ydim);
	vector<vector<Pixel>> partition_scanline(int vertex_per_raster);
	vector<vector<Pixel>> partition_with_query(int vertex_per_raster);
	QTNode *partition_qtree(const int vpr);
	QTNode *get_qtree(){
		return qtree;
	}

	void reset_grid_partition(){
		pthread_mutex_lock(&partition_lock);
		for(vector<Pixel> &rows:partitions){
			rows.clear();
		}
		partitions.clear();
		pthread_mutex_unlock(&partition_lock);
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
		assert(mbb);
		return (xval-mbb->low[0])/step_x;
	}
	inline int get_pixel_y(double yval){
		assert(mbb);
		return (yval-mbb->low[1])/step_y;
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
		area_buffer = std::max(sum,0.0);
		return area_buffer;
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

#endif /* SRC_MYPOLYGON_H_ */
