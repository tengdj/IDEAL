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

	//return box(-180,-90,180,90);
	box space;
	int nt = 2*get_num_threads()-1;
	box subspaces[nt];
#pragma omp parallel for num_threads(nt)
	for(size_t i=0;i<geometries.size();i++){
		int tid = omp_get_thread_num();
		subspaces[tid].update(*(geometries[i]->getMBB()));
	}
	for(int i=0;i<nt;i++){
		space.update(subspaces[i]);
	}
	return space;
}

class BTNode:public box{
public:
	bool isleaf(){
		return children[0]==NULL;
	}
	vector<MyPolygon *> objects;

	BTNode *children[2] = {NULL, NULL};

	BTNode(double low_x, double low_y, double high_x, double high_y){
		low[0] = low_x;
		low[1] = low_y;
		high[0] = high_x;
		high[1] = high_y;
	}
	BTNode(box m){
		low[0] = m.low[0];
		low[1] = m.low[1];
		high[0] = m.high[0];
		high[1] = m.high[1];
	}
	void split(){
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
		if(isleaf()){
			split();
		}
		for(int i=0;i<2;i++){
			children[i]->split_to(threshold);
		}
	}

	void get_leafs(vector<BTNode *> &leafs){
		if(isleaf()){
			leafs.push_back(this);
		}else{
			for(int i=0;i<2;i++){
				children[i]->get_leafs(leafs);
			}
		}
	}
	~BTNode(){
		if(!isleaf()){
			for(int i=0;i<2;i++){
				delete children[i];
			}
		}
		objects.clear();
	}
	void insert(MyPolygon *p){
		objects.push_back(p);
	}
	void insert(vector<MyPolygon *> &geometries){
		objects.insert(objects.end(), geometries.begin(), geometries.end());
	}
	size_t num_objects(){
		return objects.size();
	}

};

vector<Tile *> genschema_slc(vector<MyPolygon *> &geometries, size_t cardinality, bool data_oriented){
	vector<Tile *> schema;
	struct timeval start = get_cur_time();
	size_t num = geometries.size();
	//boost::sort::sample_sort(geometries.begin(), geometries.end(), compareCentX);
	//boost::sort::parallel_stable_sort(geometries.begin(), geometries.end(), compareCentX);
	boost::sort::block_indirect_sort(geometries.begin(), geometries.end(), compareLowX);
	//std::sort(geometries.begin(), geometries.end(), compareCentX);

	logt("sorting",start);
	size_t schema_num = std::max(num/cardinality, (size_t)1);
	schema.resize(schema_num);
	box space = profileSpace(geometries);
#pragma omp parallel for num_threads(2*get_num_threads()-1)
	for(size_t i=0;i<schema_num;i++){
		const size_t bg = i*cardinality;
		const size_t ed = std::min((i+1)*cardinality, num);
		assert(ed>bg);
		Tile *b = new Tile();
		if(data_oriented){
			for(size_t t = bg;t<ed; t++){
				b->insert(geometries[t], true);
			}
		}else{
			b->low[1] = space.low[1];
			b->high[1] = space.high[1];
			if(i==0){
				b->low[0] = space.low[0];
			}else{
				b->low[0] = geometries[bg]->getMBB()->low[0];
			}
			if(i==schema_num-1){
				b->high[0] = space.high[0];
			}else{
				b->high[0] = geometries[ed]->getMBB()->low[0];
			}
		}
		schema[i] = b;
	}
	return schema;
}


vector<Tile *> genschema_str(vector<MyPolygon *> &geometries, size_t cardinality, bool data_oriented){
	size_t part_num = geometries.size()/cardinality+1;

	size_t dimx = sqrt(part_num);
	size_t dimy = part_num/dimx;
	vector<Tile *> schema;
	struct timeval start = get_cur_time();
	size_t num = geometries.size();
	boost::sort::parallel_stable_sort(geometries.begin(), geometries.end(), compareCentX);
	schema.resize(dimx*dimy);
	box space = profileSpace(geometries);
	size_t sx = num/dimx;
	box *slice_boxes = new box[dimx];
	for(size_t x=0;x<dimx;x++){
		size_t bg_x = sx*x;
		size_t ed_x = (x+1)*sx;
		if(x==dimx){
			ed_x = num;
		}
		assert(ed_x>bg_x);
		slice_boxes[x].low[1] = space.low[1];
		slice_boxes[x].high[1] = space.high[1];
		if(x==0){
			slice_boxes[x].low[0] = space.low[0];
		}else{
			slice_boxes[x].low[0] = geometries[bg_x]->getMBB()->centroid().x;
		}
		if(x==dimx-1){
			slice_boxes[x].high[0] = space.high[0];
		}else{
			slice_boxes[x].high[0] = geometries[ed_x]->getMBB()->centroid().x;
		}
	}
#pragma omp parallel for num_threads(2*get_num_threads()-1)
	for(size_t x=0;x<dimx;x++){
		size_t bg_x = sx*x;
		size_t ed_x = std::min((x+1)*sx, num);

		//log("%d %d",bg_x,ed_x);
		boost::sort::parallel_stable_sort(geometries.begin()+bg_x, geometries.begin()+ed_x, compareCentY);

		size_t sy = (ed_x-bg_x)/dimy;
		for(size_t y=0;y<dimy;y++){
			size_t bg = y*sy+bg_x;
			size_t ed = (y+1)*sy+bg_x;
			if(y==dimy-1){
				ed = ed_x;
			}
			assert(ed>bg);
			Tile *b = new Tile();
			if(data_oriented){
				for(size_t t = bg;t<ed; t++){
					b->insert(geometries[t], true);
				}
			}else{
				b->low[0] = slice_boxes[x].low[0];
				b->high[0] = slice_boxes[x].high[0];
				if(y==0){
					b->low[1] = slice_boxes[x].low[1];
				}else{
					b->low[1] = geometries[bg]->getMBB()->centroid().y;
				}
				if(y==dimy-1){
					b->high[1] = slice_boxes[x].high[1];
				}else{
					b->high[1] = geometries[ed]->getMBB()->centroid().y;
				}
			}
			schema[y*dimx+x] = b;
		}
	}
	delete []slice_boxes;
	return schema;
}


vector<Tile *> genschema_hc(vector<MyPolygon *> &geometries, size_t cardinality, bool data_oriented){
	assert(geometries.size()>0);
	vector<Tile *> schema;
	size_t num = geometries.size();

	// make it precise enough
	size_t hcnum = geometries.size()*2;
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
	logt("calculating the hilbert curve id", start);
	if(!data_oriented){
		pthread_mutex_t hlocks[100];
		for(int i=0;i<100;i++){
			pthread_mutex_init(&hlocks[i],NULL);
		}
		uint *hccount = new uint[hcnum];
		memset((void *)hccount, 0, sizeof(uint)*hcnum);
#pragma omp parallel for num_threads(2*get_num_threads()-1)
		for(size_t i=0;i<geometries.size();i++){
			pthread_mutex_lock(&hlocks[geometries[i]->hc_id%100]);
			hccount[geometries[i]->hc_id]++;
			pthread_mutex_unlock(&hlocks[geometries[i]->hc_id%100]);
		}
		logt("assign to unit tiles", start);

		size_t thread_num = 1;//get_num_threads();
#pragma omp parallel for num_threads(thread_num)
		for(size_t i=0;i<thread_num;i++){
			size_t cur = 0;
			Tile *b = new Tile();
			uint cur_count = 0;
			size_t st = i*hcnum/thread_num;
			size_t ed = min(hcnum, (i+1)*hcnum/thread_num);
			for(size_t h=st;h<ed;h++){
				size_t x = 0;
				size_t y = 0;
				d2xy(hindex, h, x, y);
				box tb(x*space.width()/dimx+space.low[0],y*space.height()/dimy+space.low[1],(x+1)*space.width()/dimx+space.low[0],(y+1)*space.height()/dimy+space.low[1]);
				cur_count += hccount[h];
				b->update(tb);
				if(cur_count>cardinality){
#pragma omp critical
					schema.push_back(b);
					b = new Tile();
					cur_count = 0;
				}
			}
			if(cur_count>0){
#pragma omp critical
				schema.push_back(b);
			}
		}
		logt("merge adjacent tiles", start);
		delete []hccount;
	}else{
		//logt("generate hilbert curve id", start);
		boost::sort::block_indirect_sort(geometries.begin(), geometries.end(), compareHC);
		size_t schema_num = (num+cardinality-1)/cardinality;

#pragma omp parallel for num_threads(2*get_num_threads()-1)
		for(size_t i=0;i<schema_num;i++){
			size_t bg = i*cardinality;
			size_t ed = std::min((i+1)*cardinality, num);
			assert(ed>bg);

			if(bg>0){
				// skip someone that may share the same HC value with the previous tile
				size_t prev = geometries[bg-1]->hc_id;
				while(geometries[bg]->hc_id==prev){
					bg++;
				}
				if(bg>=ed){
					continue;
				}
			}

			Tile *b = new Tile();
			for(size_t t = bg;t<ed; t++){
				b->insert(geometries[t]);
			}
			// also insert some in next tile with the same
			if(ed!=geometries.size()){
				size_t cur = geometries[ed-1]->hc_id;
				while(ed<geometries.size() && geometries[ed]->hc_id==cur){
					b->insert(geometries[ed]);
					ed++;
				}
			}
#pragma omp critical
			schema.push_back(b);
		}
	}
	//logt("generate tiles", start);

	return schema;
}

vector<Tile *> genschema_fg(vector<MyPolygon *> &geometries, size_t cardinality, bool data_oriented){

	size_t part_num = geometries.size()/cardinality+1;
	size_t dimx = sqrt(part_num);
	size_t dimy = part_num/dimx;
	vector<Tile *> schema;
	box space = profileSpace(geometries);
	double sx = (space.high[0]-space.low[0])/dimx;
	double sy = (space.high[1]-space.low[1])/dimx;
	schema.resize(dimx*dimy);
	size_t cont = 0;
	for(size_t x=0;x<dimx;x++){
		for(size_t y=0;y<dimy;y++){
			Tile *t = new Tile();
			if(!data_oriented){
				t->low[0] = space.low[0]+x*sx;
				t->high[0] = space.low[0]+(x+1)*sx;
				t->low[1] = space.low[1]+y*sy;
				t->high[1] = space.low[1]+(y+1)*sy;
			}
			schema[y*dimx+x] = t;
		}
	}

	if(data_oriented){
#pragma omp parallel for num_threads(2*get_num_threads()-1)
		for(size_t i=0;i<geometries.size();i++){
			Point cent = geometries[i]->getMBB()->centroid();
			size_t x = (cent.x - space.low[0])/sx;
			size_t y = (cent.y - space.low[1])/sy;
			schema[y*dimx+x]->insert(geometries[i], true);
		}
		vector<Tile *>::iterator it = schema.begin();
		while(it!=schema.end()){
			if((*it)->objects.size()==0){
				it = schema.erase(it);
			}else{
				it++;
			}
		}
	}

	return schema;
}


vector<Tile *> genschema_qt(vector<MyPolygon *> &geometries, size_t cardinality, bool data_oriented){

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
		qtree->touch(ct, geometries[i]);
	}
	qtree->adjust(cardinality);

	vector<QTNode *> qnodes;
	qtree->get_leafs(qnodes);
	for(QTNode *qn:qnodes){
		if(data_oriented){
			if(qn->objects.size()>0){
				Tile *t = new Tile();
				for(MyPolygon *obj:qn->objects){
					t->insert(obj, true);
				}
				schema.push_back(t);
			}
		}else{
			Tile *t = new Tile(qn->mbr);
			schema.push_back(t);
		}
	}
	qnodes.clear();
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
			if(node->num_objects()>threshold){
				node->split();
				lock();
				btnodes.push(node->children[0]);
				btnodes.push(node->children[1]);
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

vector<Tile *> genschema_bsp(vector<MyPolygon *> &geometries, size_t cardinality, bool data_oriented){

	struct timeval start = get_cur_time();
	vector<Tile *> schema;
	box space = profileSpace(geometries);

	BTNode *btree = new BTNode(space);
	btree->insert(geometries);
	btnodes.push(btree);

	pthread_t threads[get_num_threads()];
	for(int i=0;i<::get_num_threads();i++){
		pthread_create(&threads[i], NULL, bsp_unit, (void *)&cardinality);
	}
	for(int i = 0; i < get_num_threads(); i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}

	vector<BTNode *> leafs;
	btree->get_leafs(leafs);

	for(BTNode *b:leafs){
		if(data_oriented){
			if(b->objects.size()>0){
				Tile *t = new Tile();
				for(void *obj:b->objects){
					t->insert((MyPolygon *)obj, true);
				}
				schema.push_back(t);
			}
		}else{
			Tile *t = new Tile(*b);
			schema.push_back(t);
		}
	}
	leafs.clear();

	delete btree;
	return schema;
}


/*
 *
 * single thread version of partitioning algorithms
 *
 * */

vector<Tile *> genschema_qt_st(vector<MyPolygon *> &geometries, size_t cardinality, bool data_oriented){

	size_t part_num = geometries.size()/cardinality+1;
	struct timeval start = get_cur_time();

	vector<Tile *> schema;
	box space = profileSpace(geometries);
	QTNode *qtree = new QTNode(space);
	size_t pnum = std::min(part_num*100, geometries.size());
	size_t max_level = (log2(pnum)/log2(4)+1);
	qtree->split_to(max_level);

	for(size_t i=0;i<geometries.size();i++){
		Point ct = geometries[i]->getMBB()->centroid();
		qtree->touch(ct, geometries[i]);
	}
	qtree->adjust(cardinality*3);

	vector<QTNode *> qnodes;
	qtree->get_leafs(qnodes);
	for(QTNode *qn:qnodes){
		if(data_oriented){
			if(qn->objects.size()>0){
				Tile *t = new Tile();
				for(MyPolygon *obj:qn->objects){
					t->insert(obj, true);
				}
				schema.push_back(t);
			}
		}else{
			Tile *t = new Tile(qn->mbr);
			schema.push_back(t);
		}
	}
	qnodes.clear();
	delete qtree;
	return schema;
}


vector<Tile *> genschema_hc_st(vector<MyPolygon *> &geometries, size_t cardinality, bool data_oriented){
	assert(geometries.size()>0);
	vector<Tile *> schema;
	size_t num = geometries.size();

	// make it precise enough
	size_t hcnum = geometries.size()*2;
	size_t hindex = log2(hcnum);
	hindex += (1+(hindex%2==0));
	hcnum = pow(2,hindex);
	size_t dimx = pow(2,hindex/2);
	size_t dimy = pow(2,hindex/2);

	box space = profileSpace(geometries);

	double sx = (space.high[0]-space.low[0])/dimx;
	double sy = (space.high[1]-space.low[1])/dimx;

	struct timeval start = get_cur_time();
	uint *hccount = new uint[hcnum];
	memset((void *)hccount, 0, sizeof(uint)*hcnum);

	for(size_t i=0;i<geometries.size();i++){
		MyPolygon *p = geometries[i];
		size_t x = (p->getMBB()->centroid().x-space.low[0])/sx;
		size_t y = (p->getMBB()->centroid().y-space.low[1])/sy;
		p->hc_id = xy2d(hindex,x,y);
		assert(p->hc_id<hcnum);
		hccount[p->hc_id]++;
	}
	logt("calculating the hilbert curve id", start);

	Tile *b = new Tile();
	uint cur_count = 0;
	for(size_t h=0;h<hcnum;h++){
		size_t x = 0;
		size_t y = 0;
		d2xy(hindex, h, x, y);
		box tb(x*space.width()/dimx+space.low[0],y*space.height()/dimy+space.low[1],(x+1)*space.width()/dimx+space.low[0],(y+1)*space.height()/dimy+space.low[1]);
		cur_count += hccount[h];
		b->update(tb);
		if(cur_count>cardinality){
			schema.push_back(b);
			b = new Tile();
			cur_count = 0;
		}
	}
	if(cur_count>0){
		schema.push_back(b);
	}
	logt("merge adjacent tiles", start);
	delete []hccount;
	return schema;
}

vector<Tile *> genschema_slc_st(vector<MyPolygon *> &geometries, size_t cardinality, bool data_oriented){
	vector<Tile *> schema;
	struct timeval start = get_cur_time();
	size_t num = geometries.size();
	size_t schema_num = std::max(num/cardinality, (size_t)1);
	size_t tiny_tile_num = 100*schema_num;
	box space = profileSpace(geometries);
	double step = space.width()/tiny_tile_num;

	vector<Tile *> tinytiles;
	for(size_t i=0;i<tiny_tile_num;i++){
		Tile *t = new Tile();
		tinytiles.push_back(t);
	}

	for(MyPolygon *p:geometries){
		size_t x = (p->getMBB()->low[0] - space.low[0])/step;
		assert(x>=0&&x<tiny_tile_num);
		tinytiles[x]->insert(p, false, false);
	}

	logt("assining",start);
	schema.resize(schema_num);

	double last_one = space.low[0];
	size_t cur_tile = 0;
	size_t last_sorted = tiny_tile_num;
	for(size_t i=0;i<schema_num;i++){
		Tile *b = new Tile();
		b->low[1] = space.low[1];
		b->high[1] = space.high[1];
		b->low[0] = last_one;
		size_t cur_count = 0;

		while(cur_tile<tiny_tile_num){
			if(cur_tile==tiny_tile_num-1 || cur_count+tinytiles[cur_tile]->objects.size()==cardinality){
				b->high[0] = space.low[0]+(cur_tile+1)*step;
				last_one = b->high[0];
				schema[i] = b;
				cur_count += tinytiles[cur_tile]->objects.size();
				cur_tile++;
				break;
			}
			if(cur_count+tinytiles[cur_tile]->objects.size()<cardinality){
				cur_count += tinytiles[cur_tile]->objects.size();
				//log("%ld %ld", cur_count,tinytiles[cur_tile]->objects.size());
				cur_tile++;
				continue;
			}
			// sort and split
			if(cur_count+tinytiles[cur_tile]->objects.size()>cardinality){
				if(last_sorted != cur_tile){
					sort(tinytiles[cur_tile]->objects.begin(),tinytiles[cur_tile]->objects.end(),compareLowX);
					last_sorted = cur_tile;
				}
				size_t fetch_count = cardinality - cur_count;
				//log("%ld %ld %ld %ld",cardinality, cur_count,tinytiles[cur_tile]->objects.size(),i);
				b->high[0] = tinytiles[cur_tile]->objects[fetch_count]->getMBB()->low[0];
				tinytiles[cur_tile]->objects.erase(tinytiles[cur_tile]->objects.begin(), tinytiles[cur_tile]->objects.begin()+fetch_count);
				last_one = b->high[0];
				schema[i] = b;
				cur_count += fetch_count;
				break;
			}
			assert(0);
		}
	}
	return schema;
}

vector<Tile *> genschema_str_st(vector<MyPolygon *> &geometries, size_t cardinality, bool data_oriented){
	const size_t ampli = 100;
	struct timeval start = get_cur_time();

	size_t part_num = geometries.size()/cardinality+1;
	const size_t dimx = sqrt(part_num);
	const size_t dimy = part_num/dimx;
	box space = profileSpace(geometries);

	vector<Tile *> schema;

	vector<vector<Tile *>> tinytiles;
	for(size_t i=0;i<dimx*ampli;i++){
		vector<Tile *> slice;
		for(size_t j=0;j<dimy*ampli;j++){
			slice.push_back(new Tile());
		}
		tinytiles.push_back(slice);
	}
	double step_x = space.width()/(dimx*ampli);
	double step_y = space.height()/(dimy*ampli);
	for(MyPolygon *p:geometries){
		size_t x = (p->getMBB()->low[0] - space.low[0])/step_x;
		size_t y = (p->getMBB()->low[1] - space.low[1])/step_y;
		assert(x>=0&&x<dimx*ampli);
		assert(y>=0&&y<dimy*ampli);
		tinytiles[x][y]->insert(p, false, false);
	}

	logt("assigning objects", start);

	space.print();

	size_t cur_slice = 0;
	double cur_right = space.low[0];
	size_t cardinality_x = geometries.size()/dimx;
	for(size_t i=0;i<dimx;i++){
		size_t cur_count = 0;
		vector<vector<Tile *>> rows;
		double cur_left = cur_right;

		while(cur_slice<dimx*ampli){

			log("%ld %ld %ld",i, cur_slice, dimx*ampli);
			// total number of objects in this slice
			size_t ct = 0;
			for(size_t j=0;j<ampli*dimy;j++){
				ct += tinytiles[cur_slice][j]->objects.size();
			}
			// accumulated slices are enough or reach the end, include the entire slice
			if(cur_slice==ampli*dimx-1 || cur_count + ct == cardinality_x){
				rows.push_back(tinytiles[cur_slice]);
				cur_right = space.low[0]+(cur_slice+1)*step_x;
				cur_slice++;
				cur_count += ct;
				//log("been here?");
				break;
			}

			// need more
			if(cur_count + ct < cardinality_x){
				//log("%ld %ld",cur_count,ct);
				rows.push_back(tinytiles[cur_slice]);
				cur_right = space.low[0]+(cur_slice+1)*step_x;
				cur_slice++;
				cur_count += ct;
				continue;
			}

			// current slice is too big, sort and split it.
			if(cur_count + ct > cardinality_x){
				//log("sort: %ld %ld",cur_count,ct);

				vector<MyPolygon *> tmp;
				vector<Tile *> left_part;
				for(size_t y=0;y<ampli*dimy;y++){
					tmp.insert(tmp.end(), tinytiles[cur_slice][y]->objects.begin(),tinytiles[cur_slice][y]->objects.end());
					tinytiles[cur_slice][y]->objects.clear();
					left_part.push_back(new Tile());
				}

				sort(tmp.begin(), tmp.end(), compareLowX);
				cur_right = tmp[cardinality_x-cur_count]->getMBB()->low[1];
				for(size_t t = 0;t<=cardinality_x-cur_count;t++){
					size_t y = (tmp[t]->getMBB()->low[1]-space.low[1])/step_y;
					left_part[y]->insert(tmp[t], false,  false);
				}
				// push back the rest to the tiny tiles
				for(size_t t = cardinality_x-cur_count+1;t<tmp.size();t++){
					size_t y = (tmp[t]->getMBB()->low[1]-space.low[1])/step_y;
					tinytiles[cur_slice][y]->insert(tmp[t], false,  false);
				}
				rows.push_back(left_part);
				break;
			}
			assert(0);
		}

		// process the other dimension with the tiles in rows
		size_t cur_row = 0;
		double last_top = space.low[1];
		const size_t card = geometries.size()/(dimx*dimy);
		for(size_t y=0;y<dimy;y++){
			Tile *b = new Tile();
			b->low[0] = cur_left;
			b->high[0] = cur_right;
			b->low[1] = last_top;
			size_t cur_count = 0;

			while(cur_row<dimy*ampli){

				size_t cty = 0;
				for(size_t x=0;x<rows.size();x++){
					cty += rows[x][cur_row]->objects.size();
				}

				if(cur_row==dimy*ampli-1 || cur_count+cty==card){
					b->high[1] = space.low[1]+(cur_row+1)*step_y;
					last_top = 	b->high[1];
					schema.push_back(b);
					cur_count += cty;
					cur_row++;
					break;
				}
				// move to next row
				if(cur_count+cty<card){
					cur_count += cty;
					//log("%ld %ld", cur_count,tinytiles[cur_tile]->objects.size());
					cur_row++;
					continue;
				}
				// sort and split
				if(cur_count+cty>card){

					vector<MyPolygon *> tmp;
					for(size_t x=0;x<rows.size();x++){
						Tile *t = rows[x][cur_row];
						tmp.insert(tmp.end(), t->objects.begin(), t->objects.end());
						t->objects.clear();
					}
					sort(tmp.begin(), tmp.end(), compareLowY);
					size_t loc = card-cur_count;
					b->high[1] = tmp[loc]->getMBB()->low[1];
					last_top = 	b->high[1];

					// push the rest into the first tile for future rounds of checking. it is unnecessary to further split them
					// into tiles.
					rows[0][cur_row]->objects.insert(rows[0][cur_row]->objects.begin(), tmp.begin()+loc+1, tmp.end());
					schema.push_back(b);
					break;
				}
				assert(0);
			}
		}

	}

	return schema;
}

class BTNode_ST:public box{
public:
	bool isleaf(){
		return children[0]==NULL;
	}
	vector<vector<Tile *>> tinytiles;
	vector<MyPolygon *> objects;

	BTNode_ST *children[2] = {NULL, NULL};

	BTNode_ST(double low_x, double low_y, double high_x, double high_y){
		low[0] = low_x;
		low[1] = low_y;
		high[0] = high_x;
		high[1] = high_y;
	}
	BTNode_ST(box m){
		low[0] = m.low[0];
		low[1] = m.low[1];
		high[0] = m.high[0];
		high[1] = m.high[1];
	}
	void split(){
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
		if(isleaf()){
			split();
		}
		for(int i=0;i<2;i++){
			children[i]->split_to(threshold);
		}
	}

	void get_leafs(vector<BTNode_ST *> &leafs){
		if(isleaf()){
			leafs.push_back(this);
		}else{
			for(int i=0;i<2;i++){
				children[i]->get_leafs(leafs);
			}
		}
	}
	~BTNode_ST(){
		if(!isleaf()){
			for(int i=0;i<2;i++){
				delete children[i];
			}
		}
		objects.clear();
	}
	void insert(MyPolygon *p){
		objects.push_back(p);
	}
	void insert(vector<MyPolygon *> &geometries){
		objects.insert(objects.end(), geometries.begin(), geometries.end());
	}
	size_t num_objects(){
		return objects.size();
	}

};

vector<Tile *> genschema_bsp_st(vector<MyPolygon *> &geometries, size_t cardinality, bool data_oriented){

	struct timeval start = get_cur_time();
	vector<Tile *> schema;
	box space = profileSpace(geometries);

	BTNode *btree = new BTNode(space);
	btree->insert(geometries, cardinality);

	vector<BTNode *> leafs;
	btree->get_leafs(leafs);

	for(BTNode *b:leafs){
		Tile *t = new Tile(*b);
		schema.push_back(t);
	}
	leafs.clear();

	delete btree;
	return schema;
}


//vector<Tile *> genschema(vector<MyPolygon *> &geometries, size_t cardinality, PARTITION_TYPE type, bool data_oriented){
//	switch(type){
//	case BSP:
//		return genschema_bsp(geometries, cardinality, data_oriented);
//	case QT:
//		return genschema_qt(geometries, cardinality, data_oriented);
//	case HC:
//		return genschema_hc(geometries, cardinality, data_oriented);
//	case SLC:
//		return genschema_slc(geometries, cardinality, data_oriented);
//	case STR:
//		return genschema_str(geometries, cardinality, data_oriented);
//	case FG:
//		return genschema_fg(geometries, cardinality, data_oriented);
//	default:
//		assert(false && "wrong partitioning type");
//	}
//	// should never reach here
//	return genschema_qt(geometries, cardinality, data_oriented);
//}

vector<Tile *> genschema(vector<MyPolygon *> &geometries, size_t cardinality, PARTITION_TYPE type, bool data_oriented){
	switch(type){
	case BSP:
		return genschema_bsp(geometries, cardinality, data_oriented);
	case QT:
		return genschema_qt_st(geometries, cardinality, data_oriented);
	case HC:
		return genschema_hc_st(geometries, cardinality, data_oriented);
	case SLC:
		return genschema_slc_st(geometries, cardinality, data_oriented);
	case STR:
		return genschema_str_st(geometries, cardinality, data_oriented);
	case FG:
		return genschema_fg(geometries, cardinality, data_oriented);
	default:
		assert(false && "wrong partitioning type");
	}
	// should never reach here
	return genschema_qt(geometries, cardinality, data_oriented);
}
