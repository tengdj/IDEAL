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
#include "util.h"
using namespace std;
const static char *multipolygon_char = "MULTIPOLYGON";
const static char *polygon_char = "POLYGON";


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

class Pixel{
public:
	int id[2];
	double low[2];
	double high[2];
	PartitionStatus status = OUT;
	PartitionStatus border[4];
	Pixel(){
		border[0] = OUT;
		border[1] = OUT;
		border[2] = OUT;
		border[3] = OUT;
	}

	vector<cross_info> crosses;
	void enter(double val, Direction d);
	void leave(double val, Direction d);
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
};

class MyPolygon{
	Pixel *mbb = NULL;
	vector<vector<Pixel>> partitions;
	double step_x = 0;
	double step_y = 0;
	bool partitioned = false;

public:
	double *boundary = NULL;
	vector<double *> internal_polygons;

	static double *read_vertices(const char *wkt, size_t &offset);
	static MyPolygon *read_polygon(const char *wkt, size_t &offset);
	static MyPolygon *gen_box(double minx,double miny,double maxx,double maxy);

	bool contain(Point &p, bool use_partition=true);
	void internal_partition();
	MyPolygon *clone();
	MyPolygon(){
	}
	~MyPolygon();
	void print_without_head();
	void print();
	int num_boundary_vertices(){
		if(!boundary){
			return 0;
		}
		return (int)(boundary[0]);
	}
	Pixel *getMBB();
	vector<vector<Pixel>> partition(int dimx, int dimy);
	inline int get_num_vertices(){
		return (int)*boundary;
	}
	inline double getx(int index){
		assert(index<*boundary);
		return boundary[index*2+1];
	}
	inline double gety(int index){
		assert(index<*boundary);
		return boundary[index*2+2];
	}
	inline int get_pixel_x(double xval){
		assert(mbb);
		return (xval-mbb->low[0])/step_x;
	}
	inline int get_pixel_y(double yval){
		assert(mbb);
		return (yval-mbb->low[1])/step_y;
	}

	static char *encode(vector<vector<Pixel>> partitions);
	static vector<vector<Pixel>> decode(char *);
	vector<Point> generate_test_points(int num);
};

class MyMultiPolygon{
	vector<MyPolygon *> polygons;
public:
	static vector<MyPolygon *> read_multipolygons(const char *wkt);
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
};

// some utility function

inline bool is_number(char ch){
	return ch=='-'||ch=='.'||(ch<='9'&&ch>='0')||ch=='e';
}

inline double read_double(const char *input, size_t &offset){
	char tmp[100];
	while(!is_number(input[offset])){
		offset++;
	}
	int index = 0;
	while(is_number(input[offset])){
		tmp[index++] = input[offset++];
	}
	tmp[index] = '\0';
	return atof(tmp);
}
inline void skip_space(const char *input, size_t &offset){
	while(input[offset]==' '||input[offset]=='\t'||input[offset]=='\n'){
		offset++;
	}
}

#endif /* SRC_MYPOLYGON_H_ */
