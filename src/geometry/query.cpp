/*
 * query.cpp
 *
 *  Created on: Jun 2, 2020
 *      Author: teng
 */


#include "MyPolygon.h"




bool MyPolygon::contain(Point &p, query_context *ctx){

	if(!mbb){
		getMBB();
	}
	if(mbb->low[0]>p.x||mbb->low[1]>p.y||
	   mbb->high[0]<p.x||mbb->high[1]<p.y){
		return false;
	}

	if(ctx&&ctx->use_partition&&partitioned){
		int part_x = this->get_pixel_x(p.x);
		int part_y = this->get_pixel_y(p.y);
		if(partitions[part_x][part_y].status==IN){
			return true;
		}
		if(partitions[part_x][part_y].status==OUT){
			return false;
		}
		bool ret = false;
		for (int i = partitions[part_x][part_y].vstart; i < partitions[part_x][part_y].vend; i++) {
			// segment i->j intersect with line y=p.y
			if ( ((gety(i)>p.y) != (gety(i+1)>p.y))){
				double a = (getx(i+1)-getx(i)) / (gety(i+1)-gety(i));
				if(p.x - getx(i) < a * (p.y-gety(i))){
					ret = !ret;
				}
			}
		}
		return ret;
	}else{
		bool ret = false;
		for (int i = 0, j = get_num_vertices()-1; i < get_num_vertices(); j = i++) {
			// segment i->j intersect with line y=p.y
			if ( ((gety(i)>p.y) != (gety(j)>p.y))){
				double a = (getx(j)-getx(i)) / (gety(j)-gety(i));
				if(p.x - getx(i) < a * (p.y-gety(i))){
					ret = !ret;
				}
			}
		}
		return ret;
	}
}


bool MyPolygon::contain_try_partition(MyPolygon *target, query_context *ctx){
	assert(partitioned&&ctx->use_partition);
	ctx->partition_determined = true;
	Pixel *a = getMBB();
	Pixel *b = target->getMBB();
	if(!a->contain(b)){
		return false;
	}
	// test all the pixels
	int txstart = this->get_pixel_x(a->low[0]);
	int txend = this->get_pixel_x(a->high[0]);
	int tystart = this->get_pixel_y(a->low[1]);
	int tyend = this->get_pixel_y(a->high[0]);
	int incount = 0;
	for(int i=txstart;i<=txend;i++){
		for(int j=tystart;j<=tyend;j++){
			if(partitions[i][j].status==OUT){
				return false;
			}
			if(partitions[i][j].status==IN){
				incount++;
			}
		}
	}
	// all is in
	if(incount==(txend-txstart+1)*(tyend-tystart+1)){
		return true;
	}

	ctx->partition_determined = false;
	return false;
}

bool MyPolygon::contain(MyPolygon *target, query_context *ctx){
	// check each vertex in target
	for(int i=0;i<target->get_num_vertices();i++){
		Point p(target->getx(i),target->gety(i));
		if(!contain(p, ctx)){
			return false;
		}
	}
	return true;
}

bool MyPolygon::contain(Pixel *target, query_context *ctx){

	Pixel *a = getMBB();
	if(!a->contain(target)){
		return false;
	}
	if(ctx&&ctx->use_partition&&partitioned){
		// test all the pixels
		int txstart = this->get_pixel_x(a->low[0]);
		int txend = this->get_pixel_x(a->high[0]);
		int tystart = this->get_pixel_y(a->low[1]);
		int tyend = this->get_pixel_y(a->high[0]);
		int incount = 0;
		for(int i=txstart;i<=txend;i++){
			for(int j=tystart;j<=tyend;j++){
				if(partitions[i][j].status==OUT){
					return false;
				}
				if(partitions[i][j].status==IN){
					incount++;
				}
			}
		}
		// all is in
		if(incount==(txend-txstart+1)*(tyend-tystart+1)){
			return true;
		}
	}

	Point p1(target->low[0], target->low[1]);
	Point p2(target->high[0], target->high[1]);
	return contain(p1, ctx)&&contain(p2, ctx);

}

bool MyPolygon::intersect(MyPolygon *target, query_context *ctx){

	Pixel *a = getMBB();
	Pixel *b = target->getMBB();
	if(!a->intersect(b)){
		return false;
	}
	if(ctx&&ctx->use_partition&&partitioned){
		// test all the pixels
		int txstart = this->get_pixel_x(a->low[0]);
		int txend = this->get_pixel_x(a->high[0]);
		int tystart = this->get_pixel_y(a->low[1]);
		int tyend = this->get_pixel_y(a->high[0]);
		int outcount = 0;
		for(int i=txstart;i<=txend;i++){
			for(int j=tystart;j<=tyend;j++){
				if(partitions[i][j].status==OUT){
					outcount++;
				}else if(partitions[i][j].status==IN){
					return true;
				}
			}
		}
		// all is out
		if(outcount==(txend-txstart+1)*(tyend-tystart+1)){
			return false;
		}
	}

	// check each vertex in target
	for(int i=0;i<target->get_num_vertices();i++){
		Point p(target->getx(i),target->gety(i));
		if(contain(p, ctx)){
			return true;
		}
	}
	return false;
}

