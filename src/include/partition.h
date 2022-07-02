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

	Tile();
	Tile(box b);
	~Tile();
	void merge(Tile *n);
	bool insert(MyPolygon *obj, bool update_box = true);
	bool insert_target(Point *obj);
	void build_index();
	vector<MyPolygon *> lookup(box *b);
	vector<MyPolygon *> lookup(Point *p);
};

vector<Tile *> genschema_str(vector<MyPolygon *> &geometries, size_t cardinality);
vector<Tile *> genschema_slc(vector<MyPolygon *> &geometries, size_t cardinality);
vector<Tile *> genschema_hc(vector<MyPolygon *> &geometries, size_t cardinality);

vector<Tile *> genschema_fg(vector<MyPolygon *> &geometries, size_t cardinality);
vector<Tile *> genschema_qt(vector<MyPolygon *> &geometries, size_t cardinality);
vector<Tile *> genschema_bsp(vector<MyPolygon *> &geometries, size_t cardinality);

vector<Tile *> genschema(vector<MyPolygon *> &geometries, size_t cardinality, PARTITION_TYPE type);
void print_tiles(vector<Tile *> &tiles);
double skewstdevratio(vector<Tile *> &tiles, int tag = 0);
PARTITION_TYPE parse_partition_type(const char *type);
bool is_data_oriented(PARTITION_TYPE);
#endif /* SRC_PARTITION_PARTITION_HPP_ */
