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

#include "MyPolygon.h"
#include "query_context.h"
#include "util.h"
#include "../index/QTree.h"
#include "../index/RTree.h"
#include "../index/hilbert_curve.hpp"
#include "partition.h"

using namespace std;

const char *partition_type_names[7] = {"str", "slc", "bos", "hc", "fg", "qt", "bsp"};


inline box profileSpace(vector<box *> &geometries){
	box space;
	for(box *p:geometries){
		space.update(*p);
	}
	return space;
}

vector<box *> genschema_fg(vector<box *> &geometries, size_t part_num){
	size_t dimx = sqrt(part_num);
	size_t dimy = part_num/dimx;
	vector<box *> schema;
	box space = profileSpace(geometries);
	double sx = (space.high[0]-space.low[0])/dimx;
	double sy = (space.high[1]-space.low[1])/dimx;

	for(size_t x=0;x<dimx;x++){
		for(size_t y=0;y<dimy;y++){
			box *p = new box();
			p->low[0] = space.low[0]+x*sx;
			p->high[0] = space.low[0]+(x+1)*sx;
			p->low[1] = space.low[1]+y*sy;
			p->high[1] = space.low[1]+(y+1)*sy;
			schema.push_back(p);
		}
	}
	return schema;
}

vector<box *> genschema_str(vector<box *> &geometries, size_t part_num){
	size_t dimx = sqrt(part_num);
	size_t dimy = part_num/dimx;
	vector<box *> schema;
	struct timeval start = get_cur_time();
	sort(geometries.begin(), geometries.end(), comparePixelX);
	logt("sort %d geometris", start, geometries.size());
	size_t num = geometries.size();

	for(size_t x=0;x<dimx;x++){
		size_t begin = (num/dimx)*x;
		size_t end = (x+1)*(num/dimx);
		if(x==dimx-1){
			end = num;
		}
		sort(geometries.begin()+begin, geometries.begin()+end, comparePixelY);
		logt("sort %d to %d (%d)", start, begin, end, end-begin);
	}


	for(size_t x=0;x<dimx;x++){
		size_t size = num/dimx;
		if(x==dimx-1){
			size = num - x*(num/dimx);
		}
		for(size_t y=0;y<dimy;y++){
			box *b = new box();
			size_t begin = x*(num/dimx)+y*(size/dimy);
			size_t end = x*(num/dimx)+(y+1)*(size/dimy);
			end = min(end,num);
			for(size_t t = begin;t<end;t++){
				b->update(*geometries[t]);
			}
			schema.push_back(b);
		}
	}
	return schema;
}

vector<box *> genschema_slc(vector<box *> &geometries, size_t part_num){
	vector<box *> schema;
	struct timeval start = get_cur_time();
	sort(geometries.begin(), geometries.end(), comparePixelX);
	logt("sort %d geometris", start, geometries.size());
	size_t num = geometries.size();

	for(size_t x=0;x<part_num;x++){
		size_t begin = x*(num/part_num);
		size_t end = (x+1)*(num/part_num);
		if(x==part_num-1){
			end = num;
		}
		box *b = new box();

		end = min(end,num);
		for(size_t t = begin;t<end;t++){
			b->update(*geometries[t]);
		}
		schema.push_back(b);
	}
	return schema;
}

typedef struct {
	box *pix;
	size_t *counter;
} intersect_counter;

bool CountIntersectNumber(box *poly, void* arg){
	intersect_counter *target = (intersect_counter *)arg;
	if(!poly->contain(*target->pix)){
		(*target->counter)++;
	}
	return true;
}

vector<box *> genschema_bos(vector<box *> &geometries, size_t part_num){
	const size_t objnum = geometries.size();
	const size_t cardinality = objnum/part_num;
	vector<box *> schema;
	box space = profileSpace(geometries);

	vector<box *> xordered;
	vector<box *> yordered;
	xordered.insert(xordered.begin(), geometries.begin(), geometries.end());
	yordered.insert(yordered.begin(), geometries.begin(), geometries.end());
	sort(xordered.begin(), xordered.end(), comparePixelX);
	sort(yordered.begin(), yordered.end(), comparePixelY);
	log("sorting");

	RTree<box *, double, 2, double> tree;
	for(box *g:geometries){
		tree.Insert(g->low, g->high, g);
	}
	log("creating index");
	size_t x_iter = 0;
	size_t y_iter = 0;

	double cur_x = xordered[x_iter]->low[0];
	double cur_y = yordered[y_iter]->low[1];

	while(x_iter<xordered.size()||y_iter<geometries.size()){
		// if cut an X slice
		size_t xsize = (x_iter+cardinality)<xordered.size()?cardinality:(xordered.size()-x_iter);
		size_t true_xsize = 0;
		box xbuffer;
		double tmp_cur_x = cur_x;
		size_t tmp_x_iter = x_iter;
		for(;tmp_x_iter<xordered.size()&&true_xsize<xsize;tmp_x_iter++){
			// already cut as a Y slice
			if(xordered[tmp_x_iter]->low[1]<cur_y){
				continue;
			}
			xbuffer.update(*xordered[tmp_x_iter]);
			tmp_cur_x = xordered[tmp_x_iter]->low[0];
			true_xsize++;
		}
		size_t xcost = 0;
		intersect_counter ct;
		ct.pix = &xbuffer;
		ct.counter = &xcost;
		tree.Search(xbuffer.low, xbuffer.high, CountIntersectNumber, (void *)&ct);

		// if cut a Y slice
		size_t ysize = (y_iter+cardinality)<yordered.size()?cardinality:(yordered.size()-y_iter);
		size_t true_ysize = 0;
		box ybuffer;
		double tmp_cur_y = cur_y;
		size_t tmp_y_iter = y_iter;
		for(;tmp_y_iter<yordered.size()&&true_ysize<ysize;tmp_y_iter++){
			// already cut as a Y slice
			if(yordered[tmp_y_iter]->low[0]<cur_x){
				continue;
			}
			ybuffer.update(*yordered[tmp_y_iter]);
			tmp_cur_y = yordered[tmp_y_iter]->low[1];
			true_ysize++;
		}
		size_t ycost = 0;
		ct.pix = &ybuffer;
		ct.counter = &ycost;
		tree.Search(ybuffer.low, ybuffer.high, CountIntersectNumber, (void *)&ct);


		if(xcost<=ycost && x_iter<xordered.size()){
			schema.push_back(new box(xbuffer));
			cur_x = tmp_cur_x;
			x_iter = tmp_x_iter;
		}else{
			schema.push_back(new box(ybuffer));
			cur_y = tmp_cur_y;
			y_iter = tmp_y_iter;
		}
		log("%ld:\t xcost %ld\t ycost %ld\t x_iter %ld\t y_iter %ld\t cur_x %.2f\t cur_y %.2f",schema.size(),xcost, ycost, x_iter, y_iter,cur_x,cur_y);

	}
	xordered.clear();
	yordered.clear();

	return schema;
}


vector<box *> genschema_hc(vector<box *> &geometries, size_t part_num){
	assert(geometries.size()>0);
	vector<box *> schema;

	size_t hcnum = part_num*1000;
	size_t hindex = log2(hcnum);
	hindex += (1+(hindex%2==0));
	hcnum = pow(2,hindex);
	size_t dimx = pow(2,hindex/2);
	size_t dimy = pow(2,hindex/2);

	size_t *cell_count = new size_t[hcnum];
	vector<box> cells;
	cells.reserve(hcnum);

	box space = profileSpace(geometries);

	double sx = (space.high[0]-space.low[0])/dimx;
	double sy = (space.high[1]-space.low[1])/dimx;

	for(box *p:geometries){
		size_t x = (p->low[0]-space.low[0])/sx;
		size_t y = (p->low[1]-space.low[1])/sy;
		size_t hc = xy2d(hindex,x,y);
		assert(hc<hcnum);
		if(cell_count[hc]==0){
			cells[hc].low[0] = p->low[0];
			cells[hc].low[1] = p->low[1];
			cells[hc].high[0] = p->high[0];
			cells[hc].high[1] = p->high[1];
		}else{
			cells[hc].update(*p);
		}
		cell_count[hc]++;
	}

	size_t avg_num = geometries.size()/part_num+1;

	box *pix = new box();
	size_t cur_num = 0;
	for(size_t i=0;i<hcnum;i++){
		if(cell_count[i]>0){
			pix->update(cells[i]);
			cur_num += cell_count[i];
		}
		if(cur_num>=avg_num){
			schema.push_back(pix);
			pix = new box();
			cur_num = 0;
		}
	}
	if(cur_num>0){
		schema.push_back(pix);
	}else{
		delete pix;
	}

	delete []cell_count;
	cells.clear();
	return schema;
}


vector<box *> genschema_qt(vector<box *> &geometries, size_t part_num){

	vector<box *> schema;
	box space = profileSpace(geometries);
	QTNode *qtree = new QTNode(space);

	size_t pnum = part_num*1000;
	size_t avg_num = geometries.size()/part_num;
	size_t max_level = (log2(pnum)/log2(4)+1);

	qtree->split_to(max_level);

	for(box *g:geometries){
		Point p(g->low[0], g->low[1]);
		qtree->touch(p);
	}
	qtree->converge(avg_num+avg_num/2);

	qtree->get_leafs(schema);

	delete qtree;
	return schema;
}

vector<box *> genschema_bsp(vector<box *> &geometries, size_t part_num){

	vector<box *> schema;
	box space = profileSpace(geometries);

	BTNode *btree = new BTNode(space);
	btree->objects.insert(btree->objects.end(), geometries.begin(), geometries.end());
	btree->split_to(geometries.size()/part_num);
	btree->get_leafs(schema);
	delete btree;
	return schema;
}

vector<box *> genschema(vector<box *> &geometries, size_t part_num, PARTITION_TYPE type){
	switch(type){
	case BSP:
		return genschema_bsp(geometries, part_num);
	case QT:
		return genschema_qt(geometries, part_num);
	case HC:
		return genschema_hc(geometries, part_num);
	case BOS:
		return genschema_bos(geometries, part_num);
	case SLC:
		return genschema_slc(geometries, part_num);
	case STR:
		return genschema_str(geometries, part_num);
	case FG:
		return genschema_fg(geometries, part_num);
	default:
		assert(false && "wrong partitioning type");
	}
	// should never reach here
	return genschema_qt(geometries, part_num);
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


