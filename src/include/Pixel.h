/*
 * Pixel.h
 *
 *  Created on: Jan 1, 2021
 *      Author: teng
 */

#ifndef SRC_GEOMETRY_PIXEL_H_
#define SRC_GEOMETRY_PIXEL_H_
#include <float.h>

#include "Point.h"
#include "geometry_computation.h"

const static char *direction_str = "lrtb";

#define BOX_GLOBAL_MIN 100000.0
#define BOX_GLOBAL_MAX -100000.0
class box{
public:
	double low[2] = {BOX_GLOBAL_MIN,BOX_GLOBAL_MIN};
	double high[2] = {BOX_GLOBAL_MAX,BOX_GLOBAL_MAX};

	box(){}

	box(box *b);
	box (double lowx, double lowy, double highx, double highy);
	bool valid();

	box get_union(box &b);
	box get_intersection(box &b);


	void update(Point &p);
	void update(box &b);

	double height();
	double width();
	double area();

	bool intersect(Point &start, Point &end);
	bool intersect(box &target);
	bool contain(box &target);
	bool contain(Point &p);

	// distance to box
	double distance(box &target, bool geography);
	double max_distance(box &target, bool geography);

	// distance to point
	double distance(Point &p, bool geography);
	double max_distance(Point &p, bool geography);

	// distance to segment
	double distance(Point &start, Point &end, bool geography);
	double max_distance(Point &start, Point &end, bool geography);

	box expand(double expand_buffer, bool geography);

	Point centroid();

	void print_vertices();
	void print();
	void to_array(Point *p);
};

/*
 *
 *	for pixel class
 *
 * */

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

class Pixels{
	uint8_t *status = nullptr;
	uint16_t *pointer = nullptr;
	pair<uint16_t, uint16_t> *edge_sequences = nullptr;   //change uint8_t to uint16_t, for many vertices case

	int len_edge_sequences = 0;
public:	
	Pixels() = default;
	Pixels(int num_vertices);
	~Pixels();

	uint16_t get_pointer(int id) {return pointer[id];}
	
	void set_status(int id, PartitionStatus status);
	PartitionStatus show_status(int id);
	void set_offset(int id, int idx){pointer[id] = idx;}
	void process_pixels_null(int x, int y);
	
	void init_edge_sequences(int num_edge_seqs);
	void add_edge(int idx, int start, int end);
	pair<uint16_t, uint16_t> get_edge_sequence(int idx){return edge_sequences[idx];}
};

/*
 * the node for RTree
 *
 * */

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
	bool validate(){
		if(is_leaf()){
			return node_element != NULL;
		}else{
			for(RTNode *ch:children){
				if(!ch->validate()){
					return false;
				}
			}
			return true;
		}
	}
};

#endif /* SRC_GEOMETRY_PIXEL_H_ */