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
	SLC = 1,
	BOS = 2,
	HC = 3,
	FG = 4,
	QT = 5,
	BSP = 6,
	PARTITION_TYPE_NUM
}PARTITION_TYPE;
extern const char *partition_type_names[7];

class Tile: public box{
	pthread_mutex_t lk;
	size_t objnum = 0;
	void lock();
	void unlock();
	static bool lookup_tree(void *, void *arg);
public:
	size_t id;
	RTree<void *, double, 2, double> tree;

	Tile();
	Tile(box b);
	~Tile();
	bool insert(box *b, void *obj);
	size_t lookup_count(box *p);
	size_t lookup_count(Point *p);
	vector<void *> lookup(box *b);
	vector<void *> lookup(Point *p);
	size_t get_objnum(){
		return objnum;
	}
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
#endif /* SRC_PARTITION_PARTITION_HPP_ */
