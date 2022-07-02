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
#include <omp.h>

#include "../index/hilbert_curve.h"
#include "MyPolygon.h"
#include "query_context.h"
#include "util.h"
#include "../index/QTree.h"
#include "../index/RTree.h"
#include "partition.h"

using namespace std;

const char *partition_type_names[7] = {"str", "slc", "hc", "fg", "qt", "bsp", "bos"};

inline bool compareLowX(MyPolygon *p1, MyPolygon *p2)
{
    return (p1->getMBB()->low[0] < p2->getMBB()->low[0]);
}

inline bool compareLowY(MyPolygon *p1, MyPolygon *p2)
{
    return (p1->getMBB()->low[1] < p2->getMBB()->low[1]);
}

inline bool compareHighX(MyPolygon *p1, MyPolygon *p2)
{
    return (p1->getMBB()->low[0] < p2->getMBB()->low[0]);
}

inline bool compareHighY(MyPolygon *p1, MyPolygon *p2)
{
    return (p1->getMBB()->low[1] < p2->getMBB()->low[1]);
}

inline bool compareCentX(MyPolygon *p1, MyPolygon *p2)
{
    return (p1->getMBB()->centroid().x < p2->getMBB()->centroid().x);
}

inline bool compareCentY(MyPolygon *p1, MyPolygon *p2)
{
    return (p1->getMBB()->centroid().y < p2->getMBB()->centroid().y);
}

inline bool compareHC(MyPolygon *p1, MyPolygon *p2)
{
    return (p1->hc_id < p2->hc_id);
}

inline box profileSpace(vector<MyPolygon *> &geometries){

	return box(-180,-90,180,90);
//	box space;
//	int nt = 2*get_num_threads()-1;
//	box subspaces[nt];
//#pragma omp parallel for num_threads(nt)
//	for(size_t i=0;i<geometries.size();i++){
//		int tid = omp_get_thread_num();
//		subspaces[tid].update(*(geometries[i]->getMBB()));
//	}
//	for(int i=0;i<nt;i++){
//		space.update(subspaces[i]);
//	}
//	return space;
}

class BTNode:public box{

public:

	bool isleaf = true;
	BTNode *children[2] = {NULL, NULL};
	vector<MyPolygon *> objects;

	BTNode(double low_x, double low_y, double high_x, double high_y){
		isleaf = true;
		low[0] = low_x;
		low[1] = low_y;
		high[0] = high_x;
		high[1] = high_y;
	}
	BTNode(box m){
		isleaf = true;
		low[0] = m.low[0];
		low[1] = m.low[1];
		high[0] = m.high[0];
		high[1] = m.high[1];
	}
	void split(){
		isleaf = false;
		// horizontally split
		if((high[0]-low[0])>(high[1]-low[1])){
			boost::sort::block_indirect_sort(objects.begin(),objects.end(),compareCentX);
			size_t half_index = objects.size()/2;
			double mid = objects[half_index]->getMBB()->low[0];
			children[0] = new BTNode(low[0], low[1], mid, high[1]);
			children[1] = new BTNode(mid, low[1], high[0],high[1]);
			children[0]->objects.insert(children[0]->objects.end(), objects.begin(), objects.begin()+half_index);
			children[1]->objects.insert(children[1]->objects.end(), objects.begin()+half_index, objects.end());
		}else{
			// vertically split
			boost::sort::block_indirect_sort(objects.begin(),objects.end(),compareCentY);
			size_t half_index = objects.size()/2;
			double mid = objects[half_index]->getMBB()->low[1];
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

vector<Tile *> genschema_slc(vector<MyPolygon *> &geometries, size_t cardinality){
	vector<Tile *> schema;
	struct timeval start = get_cur_time();
	size_t num = geometries.size();
	boost::sort::block_indirect_sort(geometries.begin(), geometries.end(), compareCentX);
	size_t schema_num = (num+cardinality-1)/cardinality;
	schema.resize(schema_num);
#pragma omp parallel for num_threads(2*get_num_threads()-1)
	for(size_t i=0;i<schema_num;i++){
		size_t bg = i*cardinality;
		size_t ed = std::min((i+1)*cardinality, num);
		assert(ed>bg);
		Tile *b = new Tile();
		for(size_t t = bg;t<ed; t++){
			b->insert(geometries[t]);
		}
		schema[i] = b;
	}
	return schema;
}


vector<Tile *> genschema_str(vector<MyPolygon *> &geometries, size_t cardinality){
	size_t part_num = geometries.size()/cardinality+1;

	size_t dimx = sqrt(part_num);
	size_t dimy = part_num/dimx;
	vector<Tile *> schema;
	struct timeval start = get_cur_time();
	size_t num = geometries.size();
	boost::sort::block_indirect_sort(geometries.begin(), geometries.end(), compareCentX);
	schema.resize(dimx*dimy);
#pragma omp parallel for num_threads(2*get_num_threads()-1)
	for(size_t x=0;x<dimx;x++){
		size_t begin = (num/dimx)*x;
		size_t end = std::min((x+1)*(num/dimx), num);
		boost::sort::block_indirect_sort(geometries.begin()+begin, geometries.begin()+end, compareCentY);

		size_t cn = end-begin;
		size_t cd = cn/dimy;
		for(size_t y=0;y<dimy;y++){
			size_t bg = y*cd;
			size_t ed = std::min((y+1)*cd, num);
			assert(ed>bg);
			Tile *b = new Tile();
			for(size_t t = bg;t<ed; t++){
				b->insert(geometries[t]);
			}
			schema[y*dimx+x] = b;
		}
	}
	return schema;
}


vector<Tile *> genschema_hc(vector<MyPolygon *> &geometries, size_t cardinality){
	assert(geometries.size()>0);
	vector<Tile *> schema;
	size_t num = geometries.size();

	// make it precise enough
	size_t hcnum = geometries.size()*10;
	size_t hindex = log2(hcnum);
	hindex += (1+(hindex%2==0));
	hcnum = pow(2,hindex);
	size_t dimx = pow(2,hindex/2);
	size_t dimy = pow(2,hindex/2);

	box space = profileSpace(geometries);

	double sx = (space.high[0]-space.low[0])/dimx;
	double sy = (space.high[1]-space.low[1])/dimx;

	struct timeval start = get_cur_time();
#pragma omp parallel for num_threads(2*get_num_threads()-1)
	for(size_t i=0;i<geometries.size();i++){
		MyPolygon *p = geometries[i];
		size_t x = (p->getMBB()->centroid().x-space.low[0])/sx;
		size_t y = (p->getMBB()->centroid().y-space.low[1])/sy;
		p->hc_id = xy2d(hindex,x,y);
		assert(p->hc_id<hcnum);
	}
	boost::sort::block_indirect_sort(geometries.begin(), geometries.end(), compareHC);
	size_t schema_num = (num+cardinality-1)/cardinality;
	schema.resize(schema_num);
#pragma omp parallel for num_threads(2*get_num_threads()-1)
	for(size_t i=0;i<schema_num;i++){
		size_t bg = i*cardinality;
		size_t ed = std::min((i+1)*cardinality, num);
		assert(ed>bg);
		Tile *b = new Tile();
		for(size_t t = bg;t<ed; t++){
			b->insert(geometries[t]);
		}
		schema[i] = b;
	}
	return schema;
}

vector<Tile *> genschema_fg(vector<MyPolygon *> &geometries, size_t cardinality){

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

vector<Tile *> genschema_qt(vector<MyPolygon *> &geometries, size_t cardinality){

	size_t part_num = geometries.size()/cardinality+1;
	struct timeval start = get_cur_time();

	vector<Tile *> schema;
	box space = profileSpace(geometries);
	QTNode *qtree = new QTNode(space);
	size_t pnum = std::min(part_num*100, geometries.size());
	size_t max_level = (log2(pnum)/log2(4)+1);
	qtree->split_to(max_level);

#pragma omp parallel for num_threads(2*get_num_threads()-1)
	for(size_t i=0;i<geometries.size();i++){
		Point ct = geometries[i]->getMBB()->centroid();
		qtree->touch(ct);
	}
	qtree->merge_objnum();
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

int active_thread = 0;
queue<BTNode *> btnodes;
void *bsp_unit(void *arg){
	size_t threshold = *((size_t *)arg);
	while(true){
		BTNode *node = NULL;
		bool has_active = false;
		lock();
		if(!btnodes.empty()){
			node = btnodes.front();
			btnodes.pop();
			active_thread++;
		}
		has_active = (active_thread>0);
		unlock();
		if(node){
			// should split
			if(node->objects.size()>threshold){
				node->split();
				lock();
				if(node->children[0]->objects.size()>threshold){
					btnodes.push(node->children[0]);
				}
				if(node->children[1]->objects.size()>threshold){
					btnodes.push(node->children[1]);
				}
				active_thread--;
				unlock();
			}
			// this thread done with this node
			lock();
			active_thread--;
			unlock();
		}else{
			// all done
			if(!has_active){
				break;
			}
			// some other threads are busy, waiting for new nodes
			usleep(1);
		}
	}
	return NULL;
}

vector<Tile *> genschema_bsp(vector<MyPolygon *> &geometries, size_t cardinality){

	struct timeval start = get_cur_time();
	vector<Tile *> schema;
	box space = profileSpace(geometries);

	BTNode *btree = new BTNode(space);
	btree->objects.insert(btree->objects.end(), geometries.begin(), geometries.end());
	btnodes.push(btree);

	pthread_t threads[get_num_threads()];
	for(int i=0;i<::get_num_threads();i++){
		pthread_create(&threads[i], NULL, bsp_unit, (void *)&cardinality);
	}
	for(int i = 0; i < get_num_threads(); i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}

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

vector<Tile *> genschema(vector<MyPolygon *> &geometries, size_t cardinality, PARTITION_TYPE type){
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

bool Tile::insert(MyPolygon *obj, bool update_box){
	lock();
	if(update_box){
		update(*obj->getMBB());
	}
	objects.push_back(obj);
	unlock();
	return true;
}

bool Tile::insert_target(Point *tgt){
	lock();
	targets.push_back(tgt);
	unlock();
	return true;
}

void Tile::build_index(){
	lock();
	for(MyPolygon *p:objects){
		tree.Insert(p->getMBB()->low, p->getMBB()->high, p);
	}
	unlock();
}

bool Tile::lookup_tree(MyPolygon *obj, void *arg){
	vector<MyPolygon *> *results = (vector<MyPolygon *> *)arg;
	results->push_back(obj);
	return true;
}

vector<MyPolygon *> Tile::lookup(box *b){
	vector<MyPolygon *> results;
	tree.Search(b->low, b->high, lookup_tree, (void *)&results);
	return results;
}

vector<MyPolygon *> Tile::lookup(Point *p){
	vector<MyPolygon *> results;
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
