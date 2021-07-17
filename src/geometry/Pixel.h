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
	int edge_id;
	cross_info(cross_type t, int e){
		type = t;
		edge_id = e;
	}
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
    double max_distance_geography(Point &p);

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

class edge_range{
public:
	int vstart = 0;
	int vend = 0;
	edge_range(int s, int e){
		vstart = s;
		vend = e;
	}
	int size(){
		return vend-vstart+1;
	}
};

class Pixel:public box{
	vector<cross_info> crosses;

public:
	unsigned short id[2];
	PartitionStatus status = OUT;
	vector<edge_range> edge_ranges;
	vector<double> intersection_nodes[4];

public:
	Pixel(){}
	void enter(double val, Direction d, int vnum);
	void leave(double val, Direction d, int vnum);
	void process_crosses(int num_edges);
	int num_edges_covered();
};


#endif /* SRC_GEOMETRY_PIXEL_H_ */
