/*
 * partition.hpp
 *
 *  Created on: Jun 8, 2022
 *      Author: teng
 */

#ifndef SRC_PARTITION_PARTITION_HPP_
#define SRC_PARTITION_PARTITION_HPP_

#include "MyPolygon.h"
#include "../index/RTree.h"

typedef enum {
	STR = 0,
	SLC,
	HC,
	FG,
	QT,
	BSP,
	PARTITION_TYPE_NUM
}PARTITION_TYPE;
extern const char *partition_type_names[7];

class Tile: public box{
	pthread_mutex_t lk;
	void lock();
	void unlock();
	static bool lookup_tree(MyPolygon *, void *arg);
public:
	size_t id;
	vector<MyPolygon *> objects;
	vector<Point *> targets;
	RTree<MyPolygon *, double, 2, double> tree;
	double indexing_latency = 0;
	double querying_latency = 0;
	Tile();
	Tile(box b);
	~Tile();
	void merge(Tile *n);
	bool insert(MyPolygon *obj, bool update_box = true, bool thread_safe = true);
	bool insert_target(Point *obj);
	void build_index();
	size_t conduct_query(bool dry_run = false);
	vector<MyPolygon *> lookup(box *b);
	vector<MyPolygon *> lookup(Point *p);
};

class Grid: public box{
public:
	size_t dimx;
	size_t dimy;
	vector<vector<Tile *>> grid;
	Grid(box space, size_t dx, size_t dy);
	~Grid();
};

vector<Tile *> genschema(vector<MyPolygon *> &geometries, size_t cardinality, PARTITION_TYPE type, bool data_oriented);
void print_tiles(vector<Tile *> &tiles);
double skewstdevratio(vector<Tile *> &tiles, int tag = 0);

inline PARTITION_TYPE parse_partition_type(const char *type){
	for(int i=0;i<7;i++){
		if(strcasecmp(type,partition_type_names[i])==0){
			return (PARTITION_TYPE)i;
		}
	}
	assert(false && "wrong partition type");
	return QT;
}

#endif /* SRC_PARTITION_PARTITION_HPP_ */
