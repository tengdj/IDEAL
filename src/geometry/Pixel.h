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


class box{
public:
	double low[2] = {100000.0,100000.0};
	double high[2] = {-100000.0,-100000.0};

	box(){}

	box(box *b){
		low[0] = b->low[0];
		high[0] = b->high[0];
		low[1] = b->low[1];
		high[1] = b->high[1];
	}

	void update(Point &p);

	double area();
	bool intersect(box &target);
	bool contain(box &target);
	bool contain(Point &p);
	double max_distance(Point &p);
	double distance(Point &p);
    double distance_geography(Point &p);

	void print_vertices();
	void print();
};

class RTNode:public box{
public:
	vector<RTNode *> children;
	void *node_element = NULL;

	RTNode(){};
	bool is_leaf(){
		return children.size()==0;
	}
	int node_count(){
		int count = 0;
		node_count(count);
		return count;
	}
	void node_count(int &count){
		count++;
		for(RTNode *px:children){
			assert(px!=this);
			px->node_count(count);
		}
	}
};

class Pixel:public box{

public:

	PartitionStatus status = OUT;
	int vstart = -1;
	int vend = -1;
	vector<cross_info> crosses[4];

public:
	Pixel(){}
//	void update(Point &p){
//		space.update(p);
//	}
//
//	double area(){
//		return space.area();
//	}
//	bool intersect(Pixel &target){
//		return space.intersect(target.space);
//	}
//	bool contain(Pixel &target){
//		return space.contain(target.space);
//	}
//	bool contain(Point &p){
//		return space.contain(p);
//	}
//	double max_distance(Point &p){
//		return space.max_distance(p);
//	}
//	double distance(Point &p){
//		return space.distance(p);
//	}
//    double distance_geography(Point &p){
//		return space.distance_geography(p);
//    }

	void enter(double val, Direction d, int vnum);
	void leave(double val, Direction d, int vnum);
};


#endif /* SRC_GEOMETRY_PIXEL_H_ */
