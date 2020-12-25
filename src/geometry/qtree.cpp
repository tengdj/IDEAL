/*
 * qtree.cpp
 *
 *  Created on: Dec 25, 2020
 *      Author: teng
 */

#include "MyPolygon.h"


QTNode::QTNode(double low_x, double low_y, double high_x, double high_y){
	mbr.low[0] = low_x;
	mbr.low[1] = low_y;
	mbr.high[0] = high_x;
	mbr.high[1] = high_y;
}
QTNode::QTNode(Pixel m){
	mbr = m;
}
void QTNode::split(){
	isleaf = false;
	double mid_x = (mbr.high[0]+mbr.low[0])/2;
	double mid_y = (mbr.high[1]+mbr.low[1])/2;
	children[bottom_left] = new QTNode(mbr.low[0],mbr.low[1],mid_x,mid_y);
	children[bottom_right] = new QTNode(mid_x,mbr.low[1],mbr.high[0],mid_y);
	children[top_left] = new QTNode(mbr.low[0],mid_y,mid_x,mbr.high[1]);
	children[top_right] = new QTNode(mid_x,mid_y,mbr.high[0],mbr.high[1]);
	for(int i=0;i<4;i++){
		children[i]->level = level+1;
	}
}
void QTNode::push(std::stack<QTNode *> &ws){
	if(!isleaf){
		for(int i=0;i<4;i++){
			ws.push(children[i]);
		}
	}
}
QTNode::~QTNode(){
	if(!isleaf){
		for(int i=0;i<4;i++){
			delete children[i];
		}
	}
}
int QTNode::size(){
	int cs = 2+4*8;
	if(!isleaf){
		for(int i=0;i<4;i++){
			cs += children[i]->size();
		}
	}
	return cs;
}
int QTNode::leaf_count(){
	if(isleaf){
		return 1;
	}else{
		int lc = 0;
		for(int i=0;i<4;i++){
			lc += children[i]->leaf_count();
		}
		return lc;
	}
}
QTNode *QTNode::retrieve(Point &p){
	assert(mbr.contain(p));
	if(this->isleaf){
		return this;
	}else{
		int offset =  2*(p.y>(mbr.low[1]+mbr.high[1])/2)+(p.x>(mbr.low[0]+mbr.high[0])/2);
		return children[offset]->retrieve(p);
	}
}

bool QTNode::within(Point &p, double within_dist){
	if(isleaf&&(interior||exterior)){
		return false;
	}
	if(mbr.distance_geography(p)<=within_dist){
		if(isleaf){
			assert(!interior&&!exterior);
			// boundary pixel
			return true;
		}else{
			// go over children
			for(int i=0;i<4;i++){
				if(children[i]->within(p, within_dist)){
					return true;
				}
			}
		}
	}
	return false;
}

void QTNode::retrieve_within(Point &p, vector<QTNode *> &candidates, double within_dist){
	if(this->mbr.distance_geography(p)<=within_dist){
		if(isleaf){
			if(!interior&&!exterior){
				candidates.push_back(this);
			}
		}else{
			for(int i=0;i<4;i++){
				children[i]->retrieve_within(p, candidates, within_dist);
			}
		}
	}
}


bool QTNode::determine_contain(Point &p){
	assert(mbr.contain(p));
	QTNode *n = retrieve(p);
	return n->interior||n->exterior;
}

bool QTNode::determine_contain(Pixel &p, bool &has_in, bool &has_ex){
	if(!mbr.intersect(&p)){
		return true;
	}
	if(isleaf){
		//border box
		if(!interior&&!exterior){
			return false;
		}
		has_ex |= exterior;
		has_in |= interior;
		//cross box
		if(has_in&&has_ex){
			return false;
		}
	}else{
		for(int i=0;i<4;i++){
			if(!children[i]->determine_contain(p, has_in, has_ex)){
				return false;
			}
		}
	}
	return true;
}

bool QTNode::determine_contain(Pixel &p){
	assert(this->level==0);
	bool has_in = false;
	bool has_out = false;
	return determine_contain(p, has_in, has_out);
}
