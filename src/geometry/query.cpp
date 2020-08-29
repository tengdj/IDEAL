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



bool MyPolygon::contain(Point p, query_context *ctx){

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
			ctx->rastor_only = true;
			return true;
		}
		if(partitions[part_x][part_y].status==OUT){
			ctx->rastor_only = true;
			return false;
		}
		ctx->rastor_only = false;
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
		ctx->edges_checked += partitions[part_x][part_y].vend-partitions[part_x][part_y].vstart;
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


bool MyPolygon::contain_try_partition(Pixel *b, query_context *ctx){
	assert(partitioned&&ctx->use_partition);
	ctx->rastor_only = true;
	Pixel *a = getMBB();
	if(!a->contain(b)){
		return false;
	}
	// test all the pixels
	int txstart = this->get_pixel_x(b->low[0]);
	int txend = this->get_pixel_x(b->high[0]);
	int tystart = this->get_pixel_y(b->low[1]);
	int tyend = this->get_pixel_y(b->high[1]);
	int incount = 0;
	int outcount = 0;
	int total = 0;

	for(int i=txstart;i<=txend;i++){
		for(int j=tystart;j<=tyend;j++){
			total++;
			if(partitions[i][j].status==OUT){
				outcount++;
			}
			if(partitions[i][j].status==IN){
				incount++;
			}
		}
	}
	// all is in
	if(incount==total){
		return true;
	}else if(outcount==total){
		return false;
	}
	//log("%d %d %d",total,incount,outcount);

	ctx->rastor_only = false;
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

bool MyPolygon::intersect_segment(Pixel *target){
	for (int i = 0; i < get_num_vertices()-1; i++) {
		// segment i->j intersect with segment
		double x1 = target->low[0];
		double x2 = target->high[0];
		double y1 = target->low[1];
		double y2 = target->high[1];
		if ( ((gety(i)>y1) != (gety(i+1)>y1))){
			double a = (getx(i+1)-getx(i)) / (gety(i+1)-gety(i));
			double b = (getx(i)*gety(i+1)-getx(i+1)*gety(i))/(gety(i+1)-gety(i));
			double xc = a*y1+b;
			if(x1<xc!=x2<xc){
				return true;
			}
		}
		if ( ((gety(i)>y2) != (gety(i+1)>y2))){
			double a = (getx(i+1)-getx(i)) / (gety(i+1)-gety(i));
			double b = (getx(i)*gety(i+1)-getx(i+1)*gety(i))/(gety(i+1)-gety(i));
			double xc = a*y2+b;
			if(x1<xc!=x2<xc){
				return true;
			}
		}
		if ( ((getx(i)>x1) != (getx(i+1)>x1))){
			double a = (gety(i+1)-gety(i)) / (getx(i+1)-getx(i));
			double b = (gety(i)*getx(i+1)-gety(i+1)*getx(i))/(getx(i+1)-getx(i));
			double yc = a*x1+b;
			if(y1<yc!=y2<yc){
				return true;
			}
		}

		if ( ((getx(i)>x2) != (getx(i+1)>x2))){
			double a = (gety(i+1)-gety(i)) / (getx(i+1)-getx(i));
			double b = (gety(i)*getx(i+1)-gety(i+1)*getx(i))/(getx(i+1)-getx(i));
			double yc = a*x2+b;
			if(y1<yc!=y2<yc){
				return true;
			}
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
		bool is_contained = false;
		if(mbr->contain(p)){
			if(partitions[part_x][part_y].status==IN){
				return 0;
			}
			if(partitions[part_x][part_y].status==BORDER){
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
			}
			is_contained = true;
		}else{
			radius = mbr->distance(p)+step_x/2;
		}


		double mindist = 1000;
		int index = 0;
		while(true){
			double sx = 360;
			double ex = -360;

			if(is_contained){
				sx = max(mbr->low[0],p.x-radius);
				ex = min(mbr->high[0],p.x+radius);
			}else{
				if(mbr->high[1]>p.y-radius&&mbr->high[1]<p.y+radius){
					double ty = sqrt(abs(radius*radius-(mbr->high[1]-p.y)*(mbr->high[1]-p.y)));
					double xval = ty+p.x;
					if(xval<mbr->high[0]&&xval>mbr->low[0]&&xval<sx){
						sx = xval;
					}
					if(xval<mbr->high[0]&&xval>mbr->low[0]&&xval>ex){
						ex = xval;
					}
					xval = p.x-ty;
					if(xval<mbr->high[0]&&xval>mbr->low[0]&&xval<sx){
						sx = xval;
					}
					if(xval<mbr->high[0]&&xval>mbr->low[0]&&xval>ex){
						ex = xval;
					}
				}
				if(mbr->low[1]>p.y-radius&&mbr->low[1]<p.y+radius){
					double ty = sqrt(abs(radius*radius-(mbr->low[1]-p.y)*(mbr->low[1]-p.y)));
					double xval = ty+p.x;
					if(xval<mbr->high[0]&&xval>mbr->low[0]&&xval<sx){
						sx = xval;
					}
					if(xval<mbr->high[0]&&xval>mbr->low[0]&&xval>ex){
						ex = xval;
					}
					xval = p.x-ty;
					if(xval<mbr->high[0]&&xval>mbr->low[0]&&xval<sx){
						sx = xval;
					}
					if(xval<mbr->high[0]&&xval>mbr->low[0]&&xval>ex){
						ex = xval;
					}
				}
				if(mbr->low[0]>p.x-radius&&mbr->low[0]<p.x+radius){
					double ty = sqrt(abs(radius*radius-(mbr->low[0]-p.x)*(mbr->low[0]-p.x)));
					double yval = ty+p.y;
					if(yval<mbr->high[1]&&yval>mbr->low[1]){
						sx = mbr->low[0];
					}
					yval = p.y-ty;
					if(yval<mbr->high[1]&&yval>mbr->low[1]){
						sx = mbr->low[0];
					}
				}
				if(mbr->high[0]>p.x-radius&&mbr->high[0]<p.x+radius){
					double ty = sqrt(abs(radius*radius-(mbr->high[0]-p.x)*(mbr->high[0]-p.x)));
					double yval = ty+p.y;
					if(yval<mbr->high[1]&&yval>mbr->low[1]){
						ex = mbr->high[0];
					}
					yval = p.y-ty;
					if(yval<mbr->high[1]&&yval>mbr->low[1]){
						ex = mbr->high[0];
					}
				}
				if(p.y>mbr->low[1]&&p.y<mbr->high[1]){
					if(p.x-radius>mbr->low[0]&&p.x-radius<sx){
						sx = p.x-radius;
					}
					if(p.x+radius<mbr->high[0]&&p.x+radius>ex){
						ex = p.x+radius;
					}
				}
			}
//			mbr->print();
//			p.print();
//			printf("%f %f %f\n",radius, sx,ex);
			if(sx<mbr->low[0]||sx>mbr->high[0]||ex<mbr->low[0]||ex>mbr->high[0]){
				break;
			}
//			assert(sx>=mbr->low[0]&&sx<=mbr->high[0]);
//			assert(ex>=mbr->low[0]&&ex<=mbr->high[0]);

//			double sy = max(mbr->low[1],p.y-radius);
//			double ey = min(mbr->high[1],p.y+radius);
//			int pixxs = get_pixel_x(sx);
//			int pixxe = get_pixel_x(ex);
//			int pixys = get_pixel_y(sy);
//			int pixye = get_pixel_y(ey);
//			for(int pixx=pixxs;pixx<=pixxe;pixx++){
//				for(int pixy=pixys;pixy<=pixye;pixy++){
//					if(partitions[pixx][pixy].status==BORDER){
//						for (int i = partitions[pixx][pixy].vstart; i <= partitions[pixx][pixy].vend; i++) {
//							double dist = point_to_segment_distance(p.x, p.y, getx(i), gety(i), getx(i+1), gety(i+1));
//							if(dist<mindist){
//								mindist = dist;
//								//cout<<"teng: "<<radius<<" "<<dist<<" "<<pixx<<" "<<pixy<<" "<<p.x<<" "<<p.y<<endl;
//							}
//						}
//					}
//				}
//			}
//			if(mindist<360){
//				return mindist;
//			}

//
//
//			if(abs(p.x-(155.091526))<0.00001 && abs(p.y-39.771422)<0.00001){
//				printf("radius: %f\n",radius);
//				printf("start x: %f\n",sx);
//				printf("end x: %f\n\n",ex);
//			}
			for(double xval=sx;xval<=ex;xval+=step_x){
				double ty = sqrt(abs(radius*radius-(xval-p.x)*(xval-p.x)));
				double yval = ty+p.y;
				if(yval>mbr->high[1]||yval<mbr->low[1]){
					yval = p.y-ty;
				}
				if(yval>mbr->high[1]||yval<mbr->low[1]){
					continue;
				}
				int pixx = get_pixel_x(xval);
				int pixy = get_pixel_y(yval);
//				printf("%f\n",radius*radius-(xval-p.x)*(xval-p.x));
//				cout<<pixx<<" "<<ty<<" "<<pixy<<" "<<partitions.size()<<" "<<partitions[0].size()<<endl;

				if(partitions[pixx][pixy].status==BORDER){
					for (int i = partitions[pixx][pixy].vstart; i <= partitions[pixx][pixy].vend; i++) {
						double dist = point_to_segment_distance(p.x, p.y, getx(i), gety(i), getx(i+1), gety(i+1));
						if(dist<mindist){
							mindist = dist;
							//cout<<"teng: "<<radius<<" "<<dist<<" "<<pixx<<" "<<pixy<<" "<<p.x<<" "<<p.y<<endl;
						}
					}
				}
			}
			if(mindist<=1000){
				break;
			}
			radius += step_x;
			if(index++>10000){
				mbr->print();
				p.print();
				printf("%f %f %f %f %ld %ld \n",radius,sx,ex,step_x,partitions.size(),partitions[0].size());
				exit(0);
			}
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

