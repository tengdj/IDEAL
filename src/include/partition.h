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
	PARTITION_TYPE_NUM,
	BOS
}PARTITION_TYPE;
extern const char *partition_type_names[7];

class Tile: public box{
	pthread_mutex_t lk;
	void lock();
	void unlock();
	static bool lookup_tree(void *, void *arg);
public:
	size_t id;
	vector<pair<box *, void *>> objects;
	vector<void *> targets;
	RTree<void *, double, 2, double> tree;

	Tile();
	Tile(box b);
	~Tile();
	void merge(Tile *n);
	bool insert(box *b, void *obj);
	bool insert_target(void *obj);
	void build_index();
	size_t lookup_count(box *p);
	size_t lookup_count(Point *p);
	vector<void *> lookup(box *b);
	vector<void *> lookup(Point *p);
};

vector<Tile *> genschema_str(vector<box *> &geometries, size_t cardinality);
vector<Tile *> genschema_slc(vector<box *> &geometries, size_t cardinality);
vector<Tile *> genschema_bos(vector<box *> &geometries, size_t cardinality);
vector<Tile *> genschema_hc(vector<box *> &geometries, size_t cardinality);

vector<Tile *> genschema_fg(vector<box *> &geometries, size_t cardinality);
vector<Tile *> genschema_qt(vector<box *> &geometries, size_t cardinality);
vector<Tile *> genschema_bsp(vector<box *> &geometries, size_t cardinality);

vector<Tile *> genschema(vector<box *> &geometries, size_t cardinality, PARTITION_TYPE type);
void print_tiles(vector<Tile *> &tiles);
double skewstdevratio(vector<Tile *> &tiles);
PARTITION_TYPE parse_partition_type(const char *type);
bool is_data_oriented(PARTITION_TYPE);
#endif /* SRC_PARTITION_PARTITION_HPP_ */
