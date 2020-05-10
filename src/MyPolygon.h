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

class MyPolygon{

public:
	double *boundary = NULL;
	vector<double *> internal_polygons;

	static double *read_vertices(const char *wkt, size_t &offset);
	static MyPolygon *read_polygon(const char *wkt, size_t &offset);
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
};

class MyMultiPolygon{
	vector<MyPolygon *> polygons;
public:
	static vector<MyPolygon *> read_multipolygons(const char *wkt);
	MyMultiPolygon(const char *wkt);
	~MyMultiPolygon();
	int num_polygons(){
		return polygons.size();
	}
	void print();
	void print_separate();
	vector<MyPolygon *> get_polygons(){
		return polygons;
	}
};



#endif /* SRC_MYPOLYGON_H_ */
