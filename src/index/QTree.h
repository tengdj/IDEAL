/*
 * QTree.h
 *
 *  Created on: Jan 1, 2021
 *      Author: teng
 */

#ifndef SRC_INDEX_QTREE_H_
#define SRC_INDEX_QTREE_H_

#include "../geometry/Pixel.h"

enum QT_Direction{
	bottom_left = 0,
	bottom_right = 1,
	top_left = 2,
	top_right = 3
};

class QTNode{

	void internal_leaf_count(int &count){
		if(isleaf){
			count++;
		}else{
			for(int i=0;i<4;i++){
				children[i]->internal_leaf_count(count);
			}
		}
	}

	void internal_border_leaf_count(int &count){
		if(isleaf){
			if(!interior&&!exterior){
				count++;
			}
		}else{
			for(int i=0;i<4;i++){
				children[i]->internal_leaf_count(count);
			}
		}
	}

public:

	bool isleaf = true;
	QTNode *children[4];
	Pixel mbr;
	bool interior = false;
	bool exterior = false;
	int level = 0;

	QTNode(double low_x, double low_y, double high_x, double high_y){
		mbr.low[0] = low_x;
		mbr.low[1] = low_y;
		mbr.high[0] = high_x;
		mbr.high[1] = high_y;
	}
	QTNode(Pixel m){
		mbr = m;
	}
	void split(){
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
	void push(std::stack<QTNode *> &ws){
		if(!isleaf){
			for(int i=0;i<4;i++){
				ws.push(children[i]);
			}
		}
	}
	~QTNode(){
		if(!isleaf){
			for(int i=0;i<4;i++){
				delete children[i];
			}
		}
	}
	int size(){
		// each nodes need 1 byte for node status
		// 2 bytes each for four children
		const int lfc = leaf_count();
		int tmp = lfc;
		int bits = 0;
		while(tmp>0){
			bits++;
			tmp /= 2;
		}
		//int bits = (bits+7)/8*8;
		bits = bits*4+2;
		return (lfc*bits+7)/8;
	}
	int leaf_count(){
		int count = 0;
		internal_leaf_count(count);
		return count;
	}
	int border_leaf_count(){
		int count = 0;
		internal_border_leaf_count(count);
		return count;
	}

	QTNode *retrieve(Point &p){
		if(!mbr.contain(p)){
			mbr.print();
			p.print();
		}
		assert(mbr.contain(p));
		if(this->isleaf){
			return this;
		}else{
			int offset =  2*(p.y>(mbr.low[1]+mbr.high[1])/2)+(p.x>(mbr.low[0]+mbr.high[0])/2);
			return children[offset]->retrieve(p);
		}
	}

	bool within(Point &p, double within_dist){
		if(isleaf&&(interior||exterior)){
			return false;
		}
		if(mbr.distance_geography(p)<=within_dist){
			if(isleaf){
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

	void retrieve_within(Point &p, vector<QTNode *> &candidates, double within_dist){
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


	bool determine_contain(Point &p){
		assert(mbr.contain(p));
		QTNode *n = retrieve(p);
		return n->interior||n->exterior;
	}

	bool determine_contain(Pixel &p, bool &has_in, bool &has_ex){
		if(!mbr.intersect(p)){
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

	bool determine_contain(Pixel &p){
		assert(this->level==0);
		bool has_in = false;
		bool has_out = false;
		return determine_contain(p, has_in, has_out);
	}


};



#endif /* SRC_INDEX_QTREE_H_ */
