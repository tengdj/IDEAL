/*
 * query.cpp
 *
 *  Created on: Jun 2, 2020
 *      Author: teng
 */


#include "MyPolygon.h"
#include <float.h>
#include <math.h>


inline double point_to_segment_distance(const double x, const double y, const double x1, double y1, const double x2, const double y2) {

  double A = x - x1;
  double B = y - y1;
  double C = x2 - x1;
  double D = y2 - y1;

  double dot = A * C + B * D;
  double len_sq = C * C + D * D;
  double param = -1;
  if (len_sq != 0) //in case of 0 length line
      param = dot / len_sq;

  double xx, yy;

  if (param < 0) {
    xx = x1;
    yy = y1;
  } else if (param > 1) {
    xx = x2;
    yy = y2;
  } else {
    xx = x1 + param * C;
    yy = y1 + param * D;
  }

  double dx = x - xx;
  double dy = y - yy;
  return sqrt(dx * dx + dy * dy);
}



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
		for (int i = 0; i < get_num_vertices()-1; i++) {
			// segment i->j intersect with line y=p.y
			if ( ((gety(i)>p.y) != (gety(i+1)>p.y))){
				double a = (getx(i+1)-getx(i)) / (gety(i+1)-gety(i));
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

double MyPolygon::distance(Point &p, query_context *ctx){
	Pixel *mbr = getMBB();
	if(ctx&&ctx->use_partition&&partitioned){
		int part_x = this->get_pixel_x(p.x);
		int part_y = this->get_pixel_y(p.y);
		double radius = min(step_x,step_y);
		if(mbr->contain(p)){
			if(partitions[part_x][part_y].status==IN){
				return 0;
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
			if(ret){
				return 0;
			}
		}else{
			radius += mbr->distance(p);
			assert(radius>0);
		}

		double sx = max(mbr->low[0],p.x-radius);
		double ex = min(mbr->high[0],p.x+radius);
		double mindist = DBL_MAX;
		while(true){
			for(double xval=sx;xval<=ex;xval+=step_x){
				double yval = sqrt(radius*radius-(xval-p.x)*(xval-p.x))+p.y;
				if(yval<=mbr->high[1]&&yval>=mbr->low[1]){
					int pixx = get_pixel_x(xval);
					int pixy = get_pixel_y(yval);
					if(partitions[pixx][pixy].status==BORDER){
						for (int i = partitions[pixx][pixy].vstart; i < partitions[pixx][pixy].vend; i++) {
							double dist = point_to_segment_distance(p.x, p.y, getx(i), gety(i), getx(i+1), gety(i+1));
							if(dist<mindist){
								mindist = dist;
							}
						}
					}
				}
			}
			if(mindist<=radius){
				break;
			}
			radius += step_x;
			//cout<<mindist<<" "<<radius<<endl;
		}
		return mindist;
	}else{
		bool ret = false;
		for (int i = 0; i < get_num_vertices()-1; i++) {
			// segment i->j intersect with line y=p.y
			if ( ((gety(i)>p.y) != (gety(i+1)>p.y))){
				double a = (getx(i+1)-getx(i)) / (gety(i+1)-gety(i));
				if(p.x - getx(i) < a * (p.y-gety(i))){
					ret = !ret;
				}
			}
		}
		if(ret){
			return 0;
		}

		double mindist = DBL_MAX;
		for (int i = 0; i < get_num_vertices()-1; i++) {
			double dist = point_to_segment_distance(p.x, p.y, getx(i), gety(i), getx(i+1), gety(i+1));
			if(dist<mindist){
				mindist = dist;
			}
		}
		return mindist;
	}
}

