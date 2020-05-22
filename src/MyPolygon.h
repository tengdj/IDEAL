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
#include "util.h"
using namespace std;
const static char *multipolygon_char = "MULTIPOLYGON";


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

class Pixel{
public:
	int id[2];
	double low[2];
	double high[2];
	PartitionStatus status = OUT;
	bool left = false;
	bool top = false;
	bool bottom = false;
	bool right = false;

	// coordinates cannot be smaller than -200
	double in_vertex = -200;
	Direction in_direction;
	bool entered = false;
	bool leaved = false;

	void enter(double val, Direction d);
	void leave(double val, Direction d);
	void process_enter_leave(Direction enter_d, double enter_val, Direction leave_d, double leave_val);
	void setin(){
		status = IN;
		left = true;
		right = true;
		top = true;
		bottom = true;
	}


};

class MyPolygon{
	MyPolygon *mbb = NULL;
public:
	double *boundary = NULL;
	vector<double *> internal_polygons;

	static double *read_vertices(const char *wkt, size_t &offset);
	static MyPolygon *read_polygon(const char *wkt, size_t &offset);
	static MyPolygon *gen_box(double minx,double miny,double maxx,double maxy);
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
	MyPolygon *getMBB();
	vector<MyPolygon *> partition(int dimx, int dimy);
	inline int get_num_vertices(){
		return (int)*boundary;
	}
	inline double get_vertex_x(int index){
		assert(index<*boundary);
		return boundary[index*2+1];
	}
	inline double get_vertex_y(int index){
		assert(index<*boundary);
		return boundary[index*2+2];
	}
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



#endif /* SRC_MYPOLYGON_H_ */
