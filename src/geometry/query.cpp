/*
 * query.cpp
 *
 *  Created on: Jun 2, 2020
 *      Author: teng
 */


#include "MyPolygon.h"
#include <float.h>
#include <math.h>


/*
 *
 * containment test
 *
 * */
bool VertexSequence::contain(Point &point) {
	bool ret = false;
	for(int i = 0,j=num_vertices-1; i < num_vertices; j=i++) {
		if( ( (y[i] >= point.y ) != (y[j] >= point.y) ) &&
				(point.x <= (x[j] - x[i]) * (point.y - y[i]) / (y[j] - y[i]) + x[i])){
			ret = !ret;
		}
	}
	return ret;
}

bool MyPolygon::contain(Point &p){

    return boundary->contain(p);
}

bool contain_rtree(RTNode *node, Point &p, query_context *ctx){
	if(!node->contain(p)){
		return false;
	}
	if(node->is_leaf()){
		Triangle *triangle = (Triangle *)node->node_element;
		struct timeval start = get_cur_time();
		bool ret = false;
		for(int i=0,j=2;i<=2;j=i++){
			Point *start = triangle->point(i);
			Point *end = triangle->point(j);
			if((start->y >= p.y ) != (end->y >= p.y)){
				double xint = (end->x - start->x) * (p.y - start->y)/ (end->y - start->y) + start->x;
				if(p.x <= xint){
					ret = !ret;
				}
			}
		}
		if(ctx->is_contain_query()){
			ctx->edge_checked.counter += 3;
			ctx->edge_checked.execution_time += get_time_elapsed(start);
		}
		return !ret;
	}
	for(RTNode *ch:node->children){
		if(contain_rtree(ch,p,ctx)){
			return true;
		}
	}
	return false;
}

bool MyPolygon::contain(Point &p, query_context *ctx){

	// the MBB may not be checked for within query
	if(!mbr->contain(p)){
		return false;
	}
	ctx->object_checked.counter++;

	struct timeval start = get_cur_time();
	if(raster){
		start = get_cur_time();
		Pixel *target = raster->get_pixel(p);
		ctx->pixel_evaluated.counter++;
		ctx->pixel_evaluated.execution_time += get_time_elapsed(start);


		if(target->status==IN){
			return true;
		}
		if(target->status==OUT){
			return false;
		}

		start = get_cur_time();
		bool ret = false;

		// checking the intersection edges in the target pixel
		for(edge_range &rg:target->edge_ranges){
			for(int i = rg.vstart; i <= rg.vend; i++) {
				int j = i+1;
				if(((boundary->y[i] >= p.y) != (boundary->y[j] >= p.y))){
					double int_x = (boundary->x[j] - boundary->x[i]) * (p.y - boundary->y[i]) / (boundary->y[j] - boundary->y[i]) + boundary->x[i];
					if(p.x <= int_x && int_x <= target->high[0]){
						ret = !ret;
					}
				}
			}
		}

		ctx->edge_checked.counter += target->num_edges_covered();
		ctx->edge_checked.execution_time += get_time_elapsed(start);

		// check the crossing nodes on the right bar
		// swap the state of ret if odd number of intersection
		// nodes encountered at the right side of the border
		struct timeval tstart = get_cur_time();
		int nc = raster->count_intersection_nodes(p);
		if(nc%2==1){
			ret = !ret;
		}
		ctx->intersection_checked.counter += nc;
		ctx->intersection_checked.execution_time += get_time_elapsed(tstart);

		ctx->border_checked.counter++;
		ctx->border_checked.execution_time += get_time_elapsed(start);
		ctx->refine_count++;

		return ret;
	}else if(qtree){
		start = get_cur_time();
		QTNode *tnode = get_qtree()->retrieve(p);
		assert(tnode->isleaf&&tnode->mbr.contain(p));
		ctx->pixel_evaluated.counter++;
		ctx->pixel_evaluated.execution_time += get_time_elapsed(start);

		//refinement step
		start = get_cur_time();
		bool contained = false;
		if(tnode->exterior){
			contained = false;
		}else if(tnode->interior){
			contained = true;
		}else{
			if(ctx->perform_refine){
				if(rtree){
					contained = contain_rtree(rtree,p,ctx);
				}else{
					contained = contain(p);
				}
			}
		}
		ctx->refine_count++;
		ctx->edge_checked.counter += get_num_vertices();
		ctx->edge_checked.execution_time += get_time_elapsed(start);
		ctx->border_checked.counter++;
		ctx->border_checked.execution_time += get_time_elapsed(start);


		// qtree return
		return contained;
	}else{

		// check the maximum enclosed rectangle (MER)
		if(mer&&mer->contain(p)){
			return true;
		}
		// check the convex hull
		if(convex_hull&&!convex_hull->contain(p)){
			return false;
		}

		// refinement step
		start = get_cur_time();
		bool contained = false;

		if(ctx->perform_refine){
			if(rtree){
				contained = contain_rtree(rtree,p,ctx);
			}else{
				contained = contain(p);
			}
		}
		ctx->refine_count++;
		ctx->edge_checked.counter += get_num_vertices();
		ctx->edge_checked.execution_time += get_time_elapsed(start);
		ctx->border_checked.counter++;
		ctx->border_checked.execution_time += get_time_elapsed(start);

		// qtree return
		return contained;
	}
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

/*
 *
 * distance calculation
 *
 * */

inline double point_to_segment_distance(const double x, const double y, const double x1, double y1, const double x2, const double y2, bool geography = false) {

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
  if(geography){
      dx = dx/degree_per_kilometer_latitude;
      dy = dy/degree_per_kilometer_longitude(y);
  }
  return sqrt(dx * dx + dy * dy);
}

double MyPolygon::distance(Point &p){
    double mindist = DBL_MAX;
    for (int i = 0; i < get_num_vertices()-1; i++) {
        double dist = point_to_segment_distance(p.x, p.y, getx(i), gety(i), getx(i+1), gety(i+1));
        if(dist<mindist){
            mindist = dist;
        }
    }
    return mindist;
}


double MyPolygon::distance_rtree(Point &p, query_context *ctx){
	assert(rtree);
	queue<RTNode *> pixq;
	for(RTNode *p:rtree->children){
		pixq.push(p);
	}
	// set this value as the MINMAXDIST
	ctx->distance = DBL_MAX;
	vector<std::pair<RTNode *, double>> tmp;
	while(pixq.size()>0){
		int s = pixq.size();
		for(int i=0;i<s;i++){
			RTNode *pix = pixq.front();
			pixq.pop();
			double mindist = pix->distance(p);
			if(mindist>ctx->distance){
				continue;
			}
			double maxdist = pix->max_distance(p);

			if(maxdist<ctx->distance){
				ctx->distance = maxdist;
			}

			tmp.push_back(std::pair<RTNode *, double>(pix, mindist));
		}

		for(std::pair<RTNode *, double> pr:tmp){
			if(pr.second<=ctx->distance){
				if(pr.first->children.size()==0){
					Triangle *triangle = (Triangle *)pr.first->node_element;
					for(int i=0,j=2;i<=2;j=i++){
						Point *start = triangle->point(i);
						Point *end = triangle->point(j);
						double dist = point_to_segment_distance(p.x, p.y, start->x, start->y, end->x, end->y);
						if(ctx->distance>dist){
							ctx->distance = dist;
						}
					}
					ctx->edge_checked.counter += 3;
				}else{
					for(RTNode *p:pr.first->children){
						pixq.push(p);
					}
				}
			}
		}

		tmp.clear();
	}

	return ctx->distance;
}

double MyPolygon::distance(Point &p, query_context *ctx){

	// distance is 0 if contained by the polygon
	query_context tmpctx = *ctx;
	tmpctx.perform_refine = true;
	bool contained = contain(p, &tmpctx);

	if(contained){
		return 0;
	}

	struct timeval start;
	ctx->object_checked.counter++;

	if(raster){
		double mbrdist = mbr->distance(p);

		//initialize the starting pixel
		Pixel *closest = raster->get_closest_pixel(p);

		int step = 0;
		double step_size = raster->get_step();

		vector<Pixel *> needprocess;
		double mindist = 10000000.0;

//		for(vector<Pixel> &rows:partitions){
//			for(Pixel &pix:rows){
//				if(pix.status==BORDER&&pix.distance_geography(p)<=ctx->distance_buffer_size){
//					ctx->refine_count++;
//					return 0;
//				}
//			}
//		}

		bool there_is_border = false;
		bool border_checked = false;
		// for filtering
		while(true){
			start = get_cur_time();
			if(step==0){
				needprocess.push_back(closest);
			}else{
				int xmin = max(0,closest->id[0]-step);
				int xmax = min(raster->get_dimx()-1,closest->id[0]+step);
				int ymin = max(0,closest->id[1]-step);
				int ymax = min(raster->get_dimy()-1,closest->id[1]+step);

				//left scan
				if(closest->id[0]-step>=0){
					for(int y=ymin;y<=ymax;y++){
						needprocess.push_back(raster->get(closest->id[0]-step,y));
					}
				}
				//right scan
				if(closest->id[0]+step<raster->get_dimx()){
					for(int y=ymin;y<=ymax;y++){
						needprocess.push_back(raster->get(closest->id[0]+step,y));
					}
				}
				//bottom scan
				if(closest->id[1]-step>=0){
					for(int x=xmin+1;x<xmax;x++){
						needprocess.push_back(raster->get(x,closest->id[1]-step));
					}
				}
				//top scan
				if(closest->id[1]+step<raster->get_dimy()){
					for(int x=xmin+1;x<xmax;x++){
						needprocess.push_back(raster->get(x,closest->id[1]+step));
					}
				}
			}
			ctx->pixel_evaluated.counter++;
			ctx->pixel_evaluated.execution_time += get_time_elapsed(start);

			// should never happen
			// all the boxes are scanned
			if(needprocess.size()==0){
				ctx->refine_count++;
				return distance(p);
			}

			for(Pixel *cur:needprocess){
				//printf("checking pixel %d %d %d\n",cur->id[0],cur->id[1],cur->status);
				if(cur->status==BORDER){
					there_is_border = true;
					start = get_cur_time();
					ctx->border_evaluated.counter++;
					// no need to check the edges of this pixel
					bool toofar = ctx->is_within_query() && cur->distance_geography(p)>ctx->distance_buffer_size;
					bool tooclose = cur->distance(p)>=mindist;
					ctx->border_evaluated.execution_time += get_time_elapsed(start);
					if(toofar||tooclose){
						continue;
					}
					border_checked = true;


					// the vector model need be checked.

					start = get_cur_time();
					ctx->border_checked.counter++;
					ctx->edge_checked.counter += cur->num_edges_covered();
					for(edge_range &rg:cur->edge_ranges){
						for (int i = rg.vstart; i <= rg.vend; i++) {
							double dist = point_to_segment_distance(p.x, p.y, getx(i), gety(i), getx(i+1), gety(i+1));
							if(dist<mindist){
								mindist = dist;
							}
						}
					}
					ctx->edge_checked.execution_time += get_time_elapsed(start);
				}
			}
			needprocess.clear();

			//printf("step:%d radius:%f mindist:%f\n",step,mbrdist+step*step_size,mindist);

			// for within query, if all firstly processed boundary pixels
			// are not within the distance, it is for sure no edge will
			// be in the specified distance
			if(ctx->is_within_query() && there_is_border != border_checked){
				return mindist;
			}

			if(mindist<mbrdist+step*step_size){
				break;
			}

			step++;
		}
		ctx->refine_count++;
		// IDEAL return
		return mindist;

	}else if(qtree){
		if(qtree->within(p, ctx->distance_buffer_size)){
			ctx->refine_count++;
			//grid indexing return
			if(ctx->perform_refine){
				if(rtree){
					return distance_rtree(p,ctx);
				}else{
					ctx->edge_checked.counter += this->get_num_vertices();
					return distance(p);
				}
			}else{
				return DBL_MAX;
			}
		}
		return DBL_MAX;
	}else{
		//checking convex
		if(ctx->is_within_query()&&convex_hull){
			double min_dist = DBL_MAX;
			for(int i=0;i<convex_hull->num_vertices-1;i++){
				double dist = ::point_to_segment_distance(p.x, p.y, convex_hull->x[i], convex_hull->y[i], convex_hull->x[i+1], convex_hull->y[i+1], true);
				if(dist<min_dist){
					min_dist = dist;
				}
			}
			if(min_dist>ctx->distance_buffer_size){
				return min_dist;
			}
		}

		ctx->refine_count++;
		//SIMPVEC return
		if(ctx->perform_refine){
			if(rtree){
				return distance_rtree(p,ctx);
			}else{
				ctx->edge_checked.counter += get_num_vertices();
				return distance(p);
			}
		}else{
			return DBL_MAX;
		}
	}
}


/*
 *
 * intersect
 *
 * */


bool MyPolygon::intersect(MyPolygon *target, query_context *ctx){

	Pixel *a = getMBB();
	Pixel *b = target->getMBB();
	if(!a->intersect(*b)){
		return false;
	}
	if(raster){
		// test all the pixels
		vector<Pixel *> covered = raster->get_pixels(b);
		int outcount = 0;
		int incount = 0;
		for(Pixel *pix:covered){
			if(pix->status==OUT){
				outcount++;
			}else if(pix->status==IN){
				incount++;
			}
		}
		int total = covered.size();
		covered.clear();
		// all is out
		if(outcount==total){
			return false;
		}
		if(incount==total){
			return true;
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




