/*
 * partition.cpp
 *
 *  Created on: May 28, 2022
 *      Author: teng
 *
 *
 *  this tool helps to partition the space into N*M tiles based on the distribution of the points
 *
 */

#include <queue>
#include <fstream>
#include <cmath>
#include <vector>
#include <boost/sort/sort.hpp>

#include "../index/hilbert_curve.h"
#include "MyPolygon.h"
#include "query_context.h"
#include "util.h"
#include "../index/QTree.h"
#include "../index/RTree.h"
#include "partition.h"

using namespace std;

const char *partition_type_names[7] = {"str", "slc", "hc", "fg", "qt", "bsp", "bos"};

inline bool comparePixelLowX(box *p1, box *p2)
{
    return (p1->low[0] < p2->low[0]);
}

inline bool comparePixelLowY(box *p1, box *p2)
{
    return (p1->low[1] < p2->low[1]);
}

inline bool comparePixelHighX(box *p1, box *p2)
{
    return (p1->low[0] < p2->low[0]);
}

inline bool comparePixelHighY(box *p1, box *p2)
{
    return (p1->low[1] < p2->low[1]);
}

inline bool comparePixelCentX(box *p1, box *p2)
{
    return (p1->centroid().x < p2->centroid().x);
}

inline bool comparePixelCentY(box *p1, box *p2)
{
    return (p1->centroid().y < p2->centroid().y);
}


inline box profileSpace(vector<box *> &geometries){
	box space;
	for(box *p:geometries){
		space.update(*p);
	}
	return space;
}

class BTNode:public box{

public:

	bool isleaf = true;
	BTNode *children[2];
	vector<box *> objects;

	BTNode(double low_x, double low_y, double high_x, double high_y){
		isleaf = true;
		low[0] = low_x;
		low[1] = low_y;
		high[0] = high_x;
		high[1] = high_y;
	}
	BTNode(box m){
		isleaf = true;
		*this = m;
	}
	void split(){
		isleaf = false;
		// horizontally split
		if((high[0]-low[0])>(high[1]-low[1])){
			boost::sort::block_indirect_sort(objects.begin(),objects.end(),comparePixelLowX);
			size_t half_index = objects.size()/2;
			double mid = objects[half_index]->low[0];
			children[0] = new BTNode(low[0], low[1], mid, high[1]);
			children[1] = new BTNode(mid, low[1], high[0],high[1]);
			children[0]->objects.insert(children[0]->objects.end(), objects.begin(), objects.begin()+half_index);
			children[1]->objects.insert(children[1]->objects.end(), objects.begin()+half_index, objects.end());
		}else{
			// vertically split
			boost::sort::block_indirect_sort(objects.begin(),objects.end(),comparePixelLowY);
			size_t half_index = objects.size()/2;
			double mid = objects[half_index]->low[1];
			children[0] = new BTNode(low[0], low[1], high[0], mid);
			children[1] = new BTNode(low[0], mid, high[0],high[1]);
			children[0]->objects.insert(children[0]->objects.end(), objects.begin(), objects.begin()+half_index);
			children[1]->objects.insert(children[1]->objects.end(), objects.begin()+half_index, objects.end());
		}
	}
	void split_to(const size_t threshold){

		if(objects.size()<=threshold){
			return;
		}
		if(isleaf){
			split();
		}
		for(int i=0;i<2;i++){
			children[i]->split_to(threshold);
		}
	}

	void get_leafs(vector<box *> &leafs){
		if(isleaf){
			leafs.push_back(new box(*this));
		}else{
			for(int i=0;i<2;i++){
				children[i]->get_leafs(leafs);
			}
		}
	}
	~BTNode(){
		if(!isleaf){
			for(int i=0;i<2;i++){
				delete children[i];
			}
		}
		objects.clear();
	}
};

vector<Tile *> genschema_str(vector<box *> &geometries, size_t cardinality){
	size_t part_num = geometries.size()/cardinality+1;

	size_t dimx = sqrt(part_num);
	size_t dimy = part_num/dimx;
	vector<Tile *> schema;
	struct timeval start = get_cur_time();
	size_t num = geometries.size();
	boost::sort::block_indirect_sort(geometries.begin(), geometries.end(), comparePixelLowX);

	for(size_t x=0;x<dimx;x++){
		size_t begin = (num/dimx)*x;
		size_t end = (x+1)*(num/dimx);
		if(x==dimx-1){
			end = num;
		}
		boost::sort::block_indirect_sort(geometries.begin()+begin, geometries.begin()+end, comparePixelLowY);

		size_t cur = begin;
		while(cur<end){
			Tile *b = new Tile();
			for(size_t t = 0;t<cardinality && cur<end; t++, cur++){
				box *obj = geometries[cur];
				b->insert(obj, (void *)obj);
			}
			schema.push_back(b);
		}
	}
	return schema;
}

vector<Tile *> genschema_slc(vector<box *> &geometries, size_t cardinality){

	vector<Tile *> schema;
	struct timeval start = get_cur_time();
	size_t num = geometries.size();
	boost::sort::block_indirect_sort(geometries.begin(), geometries.end(), comparePixelLowX);

	size_t cur = 0;
	while(cur<num){
		Tile *b = new Tile();
		for(size_t t = 0;t<cardinality && cur<num; t++, cur++){
			box *obj = geometries[cur];
			b->insert(obj, (void *)obj);
		}
		schema.push_back(b);
	}
	return schema;
}

vector<Tile *> genschema_hc(vector<box *> &geometries, size_t cardinality){
	assert(geometries.size()>0);

	size_t part_num = geometries.size()/cardinality+1;

	vector<Tile *> schema;

	size_t hcnum = geometries.size()/10;
	size_t hindex = log2(hcnum);
	hindex += (1+(hindex%2==0));
	hcnum = pow(2,hindex);
	size_t dimx = pow(2,hindex/2);
	size_t dimy = pow(2,hindex/2);

	vector<Tile *> cells;
	cells.reserve(hcnum);
	for(size_t i=0;i<hcnum;i++){
		cells[i] = new Tile();
	}

	box space = profileSpace(geometries);

	double sx = (space.high[0]-space.low[0])/dimx;
	double sy = (space.high[1]-space.low[1])/dimx;

	for(box *p:geometries){
		size_t x = (p->low[0]-space.low[0])/sx;
		size_t y = (p->low[1]-space.low[1])/sy;
		size_t hc = xy2d(hindex,x,y);
		assert(hc<hcnum);
		cells[hc]->insert(p, (void *)p);
	}

	Tile *pix = new Tile();
	for(size_t i=0;i<hcnum;i++){
		pix->merge(cells[i]);
		if(pix->objects.size()>=cardinality){
			schema.push_back(pix);
			pix = new Tile();
		}
		delete cells[i];
	}
	cells.clear();

	if(pix->objects.size()>0){
		schema.push_back(pix);
	}else{
		delete pix;
	}

	return schema;
}

vector<Tile *> genschema_fg(vector<box *> &geometries, size_t cardinality){

	size_t part_num = geometries.size()/cardinality+1;
	size_t dimx = sqrt(part_num);
	size_t dimy = part_num/dimx;
	vector<Tile *> schema;
	box space = profileSpace(geometries);
	double sx = (space.high[0]-space.low[0])/dimx;
	double sy = (space.high[1]-space.low[1])/dimx;

	for(size_t x=0;x<dimx;x++){
		for(size_t y=0;y<dimy;y++){
			Tile *t = new Tile();
			t->low[0] = space.low[0]+x*sx;
			t->high[0] = space.low[0]+(x+1)*sx;
			t->low[1] = space.low[1]+y*sy;
			t->high[1] = space.low[1]+(y+1)*sy;
			schema.push_back(t);
		}
	}

	return schema;
}

vector<Tile *> genschema_qt(vector<box *> &geometries, size_t cardinality){

	size_t part_num = geometries.size()/cardinality+1;

	vector<Tile *> schema;
	box space = profileSpace(geometries);
	QTNode *qtree = new QTNode(space);

	size_t pnum = std::min(part_num*100, geometries.size());
	size_t max_level = (log2(pnum)/log2(4)+1);

	qtree->split_to(max_level);
	for(box *g:geometries){
		Point ct = g->centroid();
		qtree->touch(ct);
	}
	qtree->converge(3*cardinality);

	vector<box *> leafs;
	qtree->get_leafs(leafs);

	for(box *b:leafs){
		schema.push_back(new Tile(*b));
		delete b;
	}
	leafs.clear();

	delete qtree;
	return schema;
}

vector<Tile *> genschema_bsp(vector<box *> &geometries, size_t cardinality){

	vector<Tile *> schema;
	box space = profileSpace(geometries);

	BTNode *btree = new BTNode(space);
	btree->objects.insert(btree->objects.end(), geometries.begin(), geometries.end());
	btree->split_to(cardinality);

	vector<box *> leafs;
	btree->get_leafs(leafs);

	for(box *b:leafs){
		schema.push_back(new Tile(*b));
		delete b;
	}
	leafs.clear();

	delete btree;
	return schema;
}

vector<Tile *> genschema(vector<box *> &geometries, size_t cardinality, PARTITION_TYPE type){
	switch(type){
	case BSP:
		return genschema_bsp(geometries, cardinality);
	case QT:
		return genschema_qt(geometries, cardinality);
	case HC:
		return genschema_hc(geometries, cardinality);
	case SLC:
		return genschema_slc(geometries, cardinality);
	case STR:
		return genschema_str(geometries, cardinality);
	case FG:
		return genschema_fg(geometries, cardinality);
	default:
		assert(false && "wrong partitioning type");
	}
	// should never reach here
	return genschema_qt(geometries, cardinality);
}

PARTITION_TYPE parse_partition_type(const char *type){
	for(int i=0;i<7;i++){
		if(strcasecmp(type,partition_type_names[i])==0){
			return (PARTITION_TYPE)i;
		}
	}
	assert(false && "wrong partition type");
	return QT;
}

bool is_data_oriented(PARTITION_TYPE pt){
	switch(pt){
	case STR:
	case SLC:
	case HC:
		return true;
	case FG:
	case QT:
	case BSP:
		return false;
	default:
		assert(false && "wrong partitioning type");
		return false;
	}
}

Tile::Tile(){
	pthread_mutex_init(&lk, NULL);
	id = 0;
}

Tile::Tile(box b){
	pthread_mutex_init(&lk, NULL);
	id = 0;
	low[0] = b.low[0];
	low[1] = b.low[1];
	high[0] = b.high[0];
	high[1] = b.high[1];

}

Tile::~Tile(){
	targets.clear();
	objects.clear();
}

void Tile::lock(){
	pthread_mutex_lock(&lk);
}

void Tile::unlock(){
	pthread_mutex_unlock(&lk);
}

void Tile::merge(Tile *n){
	if(n->objects.size()==0&&n->targets.size()==0){
		return;
	}
	lock();
	update(*n);
	objects.insert(objects.end(),n->objects.begin(), n->objects.end());
	targets.insert(targets.end(),n->targets.begin(), n->targets.end());
	unlock();
}

bool Tile::insert(box *b, void *obj, bool update_box){
	lock();
	if(update_box){
		update(*b);
	}
	objects.push_back(pair<box *, void *>(b, obj));
	unlock();
	return true;
}

bool Tile::insert_target(void *tgt){
	lock();
	targets.push_back(tgt);
	unlock();
	return true;
}

void Tile::build_index(){
	lock();
	for(pair<box *, void *> &p:objects){
		tree.Insert(p.first->low, p.first->high, p.second);
	}
	unlock();
}

bool Tile::lookup_tree(void *obj, void *arg){
	vector<void *> *results = (vector<void *> *)arg;
	results->push_back(obj);
	return true;
}

vector<void *> Tile::lookup(box *b){
	vector<void *> results;
	tree.Search(b->low, b->high, lookup_tree, (void *)&results);
	return results;
}

vector<void *> Tile::lookup(Point *p){
	vector<void *> results;
	tree.Search((double *)p, (double *)p, lookup_tree, (void *)&results);
	return results;
}

void print_tiles(vector<Tile *> &boxes){
	MyMultiPolygon *cboxes = new MyMultiPolygon();
	for(Tile *p:boxes){
		MyPolygon *m = MyPolygon::gen_box(*p);
		cboxes->insert_polygon(m);
	}
	cboxes->print();
	delete cboxes;
}

double skewstdevratio(vector<Tile *> &tiles, int tag){
	if(tiles.size()==0){
		return 0;
	}
	size_t total = 0;
	for(Tile *t:tiles){
		size_t obj_num = 0;
		switch(tag){
		case 0:
			obj_num += t->objects.size();
			break;
		case 1:
			obj_num += t->targets.size();
			break;
		case 2:
			obj_num += t->objects.size();
			obj_num += t->targets.size();
			break;
		default:
			assert(false&&"can only be 0 1 2");
		}
		total += obj_num;
	}
	double avg = 1.0*total/tiles.size();
	double st = 0.0;
	for(Tile *t:tiles){
		size_t obj_num = 0;
		switch(tag){
		case 0:
			obj_num += t->objects.size();
			break;
		case 1:
			obj_num += t->targets.size();
			break;
		case 2:
			obj_num += t->objects.size();
			obj_num += t->targets.size();
			break;
		default:
			assert(false&&"can only be 0 1 2");
		}
		st += (obj_num-avg)*(obj_num-avg)/tiles.size();
	}
	return sqrt(st)/avg;
}
