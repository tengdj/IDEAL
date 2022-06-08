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

using namespace std;

Pixel profileSpace(vector<Pixel *> &geometries){

	Pixel space;
	for(Pixel *p:geometries){
		space.update(*p);
	}
	return space;
}

vector<Pixel *> genschema_fg(vector<Pixel *> &geometries, size_t part_num){
	size_t dimx = sqrt(part_num);
	size_t dimy = part_num/dimx;
	vector<Pixel *> schema;
	Pixel space = profileSpace(geometries);
	double sx = (space.high[0]-space.low[0])/dimx;
	double sy = (space.high[1]-space.low[1])/dimx;

	for(size_t x=0;x<dimx;x++){
		for(size_t y=0;y<dimy;y++){
			Pixel *p = new Pixel();
			p->low[0] = space.low[0]+x*sx;
			p->high[0] = space.low[0]+(x+1)*sx;
			p->low[1] = space.low[1]+y*sy;
			p->high[1] = space.low[1]+(y+1)*sy;
			schema.push_back(p);
		}
	}
	return schema;
}

vector<Pixel *> genschema_str(vector<Pixel *> &geometries, size_t part_num){
	size_t dimx = sqrt(part_num);
	size_t dimy = part_num/dimx;
	vector<Pixel *> schema;
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
			Pixel *b = new Pixel();
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

vector<Pixel *> genschema_slc(vector<Pixel *> &geometries, size_t part_num){
	vector<Pixel *> schema;
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
		Pixel *b = new Pixel();

		end = min(end,num);
		for(size_t t = begin;t<end;t++){
			b->update(*geometries[t]);
		}
		schema.push_back(b);
	}
	return schema;
}


vector<Pixel *> genschema_hc(vector<Pixel *> &geometries, size_t part_num){
	assert(geometries.size()>0);
	vector<Pixel *> schema;

	size_t hcnum = part_num*1000;
	size_t hindex = log2(hcnum);
	hindex += (1+(hindex%2==0));
	hcnum = pow(2,hindex);
	size_t dimx = pow(2,hindex/2);
	size_t dimy = pow(2,hindex/2);

	size_t *cell_count = new size_t[hcnum];
	vector<Pixel> cells;
	cells.reserve(hcnum);

	Pixel space = profileSpace(geometries);

	double sx = (space.high[0]-space.low[0])/dimx;
	double sy = (space.high[1]-space.low[1])/dimx;

	for(Pixel *p:geometries){
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

	Pixel *pix = new Pixel();
	size_t cur_num = 0;
	for(size_t i=0;i<hcnum;i++){
		if(cell_count[i]>0){
			pix->update(cells[i]);
			cur_num += cell_count[i];
		}
		if(cur_num>=avg_num){
			schema.push_back(pix);
			pix = new Pixel();
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


vector<Pixel *> genschema_qt(vector<Pixel *> &geometries, size_t part_num){

	vector<Pixel *> schema;
	Pixel space = profileSpace(geometries);
	QTNode *qtree = new QTNode(space);

	size_t pnum = part_num*1000;
	size_t avg_num = geometries.size()/part_num;
	size_t max_level = (log2(pnum)/log2(4)+1);

	qtree->split_to(max_level);

	for(Pixel *g:geometries){
		Point p(g->low[0], g->low[1]);
		qtree->touch(p);
	}
	qtree->converge(avg_num+avg_num/2);

	qtree->get_leafs(schema);

	delete qtree;
	return schema;
}




vector<Pixel *> genschema_bsp(vector<Pixel *> &geometries, size_t part_num){

	vector<Pixel *> schema;
	Pixel space = profileSpace(geometries);

	BTNode *btree = new BTNode(space);
	btree->objects.insert(btree->objects.end(), geometries.begin(), geometries.end());
	btree->split_to(geometries.size()/part_num);
	btree->get_leafs(schema);
	delete btree;
	return schema;
}




