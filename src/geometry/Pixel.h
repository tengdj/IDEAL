/*
 * Pixel.h
 *
 *  Created on: Jan 1, 2021
 *      Author: teng
 */

#ifndef SRC_GEOMETRY_PIXEL_H_
#define SRC_GEOMETRY_PIXEL_H_

#include "Point.h"
#include <float.h>

const static char *direction_str = "lrtb";

enum PartitionStatus{
	OUT = 0,
	BORDER = 1,
	IN = 2
};

enum Direction{
	LEFT = 0,
	RIGHT = 1,
	TOP = 2,
	BOTTOM
};

enum cross_type{
	ENTER = 0,
	LEAVE = 1
};

class cross_info{
public:
	cross_type type;
	double vertex;
};

class range{
public:
	double m_min = DBL_MAX;
	double m_max = 0;
};

class Pixel{
public:
	int id[2];
	vector<Pixel *> children;
	void *node_element = NULL;
	double low[2] = {100000.0,100000.0};
	double high[2] = {-100000.0,-100000.0};
	PartitionStatus status = OUT;

	int vstart = -1;
	int vend = -1;
	Pixel(){

	}

	~Pixel(){
		for(Pixel *p:children){
			delete p;
		}
	}

	Pixel(Pixel *p){
		low[0] = p->low[0];
		low[1] = p->low[1];
		high[0] = p->high[0];
		high[1] = p->high[1];
	}

	void update(Point &p);
	bool intersect(Pixel *target);
	bool contain(Pixel *target);
	bool contain(Point &p);
	double max_distance(Point &p);
	double area();
	double distance(Point &p);
	range distance_range(Point &p);
    double distance_geography(Point &p);

	vector<cross_info> crosses[4];
	void enter(double val, Direction d, int vnum);
	void leave(double val, Direction d, int vnum);

	void print();
	void print_polygon();
	void print_vertices();

};


#endif /* SRC_GEOMETRY_PIXEL_H_ */
