/*
 * query.cpp
 *
 *  Created on: Jun 2, 2020
 *      Author: teng
 */


#include "MyPolygon.h"
#include <float.h>
#include <math.h>
#include "geometry_computation.h"

/*
 *
 * containment test
 *
 * */
bool VertexSequence::contain(Point &point) {
	bool ret = false;
	for(int i = 0,j=num_vertices-1; i < num_vertices; j=i++) {
		if( ( (p[i].y >= point.y ) != (p[j].y >= point.y) ) &&
				(point.x <= (p[j].x - p[i].x) * (point.y - p[i].y) / (p[j].y - p[i].y) + p[i].x)){
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
		uint edge_count = 0;
		for(edge_range &rg:target->edge_ranges){
			for(int i = rg.vstart; i <= rg.vend; i++) {
				int j = i+1;
				if(((boundary->p[i].y >= p.y) != (boundary->p[j].y >= p.y))){
					double int_x = (boundary->p[j].x - boundary->p[i].x) * (p.y - boundary->p[i].y) / (boundary->p[j].y - boundary->p[i].y) + boundary->p[i].x;
					if(p.x <= int_x && int_x <= target->high[0]){
						ret = !ret;
					}
				}
			}
			edge_count += rg.size();
		}
		ctx->edge_checked.counter += edge_count;
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

		return contained;
	}
}

bool MyPolygon::contain(MyPolygon *target, query_context *ctx){
	if(!getMBB()->contain(*target->getMBB())){
		//log("mbb do not contain");
		return false;
	}

	if(ctx->use_grid && raster){
		vector<Pixel *> pxs = raster->retrieve_pixels(target->getMBB());
		int etn = 0;
		int itn = 0;
		vector<Pixel *> bpxs;
		for(Pixel *p:pxs){
			if(p->status==OUT){
				etn++;
			}else if(p->status==IN){
				itn++;
			}else{
				bpxs.push_back(p);
			}
		}
		//log("%d %d %d",etn,itn,pxs.size());

		if(etn == pxs.size()){
			return false;
		}
		if(itn == pxs.size()){
			return true;
		}

		struct timeval start = get_cur_time();
		if(target->raster){
			vector<Pixel *> bpxs2;
			for(Pixel *p:bpxs){
				bpxs2 = target->raster->retrieve_pixels(p);
				for(Pixel *p2:bpxs2){
					for(edge_range &r:p->edge_ranges){
						for(edge_range &r2:p2->edge_ranges){
							if(segment_intersect_batch(boundary->p+r.vstart, target->boundary->p+r2.vstart, r.size(), r2.size())){
								return false;
							}
						}
					}
				}
				bpxs2.clear();
			}
		}else{
			for(Pixel *p:bpxs){
				for(edge_range &r:p->edge_ranges){
					if(segment_intersect_batch(this->boundary->p+r.vstart, target->boundary->p, r.size(), target->boundary->num_vertices)){
						logt("%ld boundary %d(%ld) %d(%ld)",start,bpxs.size(),getid(),this->get_num_vertices(),target->getid(), target->get_num_vertices());
						return false;
					}
				}
			}
		}
		bpxs.clear();
	} else {
		// checking the bounding box of the target first
		Point mbb_vertices[5];
		target->mbr->to_array(mbb_vertices);
		// no intersection between this polygon and the mbr of the target polygon
		if(!segment_intersect_batch(boundary->p, mbb_vertices, boundary->num_vertices, 5)){
			// the target must be contained as its mbr is contained
			if(contain(mbb_vertices[0], ctx)){
				return true;
			}
		}
		// otherwise, checking all the edges to make sure no intersection
		if(segment_intersect_batch(boundary->p, target->boundary->p, boundary->num_vertices, target->boundary->num_vertices)){
			return false;
		}
	}
	//log("do not intersect");
	Point p(target->getx(0),target->gety(0));
	return contain(p, ctx);

}

/*
 * distance related
 *
 * */

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
			double mindist = pix->distance(p,ctx->geography);
			if(mindist>ctx->distance){
				continue;
			}
			double maxdist = pix->max_distance(p,ctx->geography);

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
						double dist = point_to_segment_distance(p, *triangle->point(i), *triangle->point(j),ctx->geography);
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
	double mindist = getMBB()->max_distance(p, ctx->geography);
	query_context tmpctx = *ctx;
	tmpctx.perform_refine = true;
	struct timeval query_start = get_cur_time();
	bool contained = contain(p, &tmpctx);
	ctx->contain_check.execution_time += get_time_elapsed(query_start,true);
	ctx->contain_check.counter++;
	if(contained){
		return 0;
	}

	ctx->object_checked.counter++;

	if(raster){
		double mbrdist = mbr->distance(p,ctx->geography);

		//initialize the starting pixel
		Pixel *closest = raster->get_closest_pixel(p);
		int step = 0;
		double step_size = raster->get_step(ctx->geography);
		vector<Pixel *> needprocess;

		bool there_is_border = false;
		bool border_checked = false;
		while(true){
			struct timeval start = get_cur_time();
			if(step==0){
				needprocess.push_back(closest);
			}else{
				needprocess = raster->expand_radius(closest, step);
			}
			// should never happen
			// all the boxes are scanned
			if(needprocess.size()==0){
				//assert(false&&"should not evaluated all boxes");
				ctx->refine_count++;
				return point_to_segment_distance_batch(p, boundary->p, boundary->num_vertices,ctx->geography);
			}
			ctx->pixel_evaluated.counter += needprocess.size();
			ctx->pixel_evaluated.execution_time += get_time_elapsed(start);

			for(Pixel *cur:needprocess){
				//printf("checking pixel %d %d %d\n",cur->id[0],cur->id[1],cur->status);
				if(cur->is_boundary()){
					there_is_border = true;
					start = get_cur_time();
					ctx->border_evaluated.counter++;
					// no need to check the edges of this pixel
					bool toofar = ctx->is_within_query() && cur->distance(p, ctx->geography)>ctx->distance_buffer_size;
					bool tooclose = cur->distance(p,ctx->geography)>=mindist;
					ctx->border_evaluated.execution_time += get_time_elapsed(start);
					if(toofar||tooclose){
						continue;
					}
					border_checked = true;

					// the vector model need be checked.
					start = get_cur_time();
					ctx->border_checked.counter++;
					int edge_num = 0;
					for(edge_range &rg:cur->edge_ranges){
						for (int i = rg.vstart; i <= rg.vend; i++) {
							double dist = point_to_segment_distance(p, *get_point(i), *get_point(i+1),ctx->geography);
							if(dist<mindist){
								mindist = dist;
							}
						}
						edge_num += rg.size();
					}
					ctx->edge_checked.counter += cur->num_edges_covered();
					ctx->edge_checked.execution_time += get_time_elapsed(start);
				}
			}
			//printf("point to polygon distance - step:%d #pixels:%ld radius:%f mindist:%f\n",step,needprocess.size(),mbrdist+step*step_size,mindist);
			needprocess.clear();
			// for within query, if all firstly processed boundary pixels
			// are not within the distance, it is for sure no edge will
			// be in the specified distance
			if(ctx->is_within_query() && there_is_border != border_checked){
				return mindist;
			}

			if(mindist<=mbrdist+step*step_size){
				break;
			}

			step++;
		}
		ctx->refine_count++;
		// IDEAL return
		return mindist;

	}else if(qtree){
		if(ctx->is_within_query()&&
				!qtree->within(p, ctx->distance_buffer_size)){
			return DBL_MAX;
		}

		ctx->refine_count++;
		if(ctx->perform_refine){
			ctx->edge_checked.counter += this->get_num_vertices();
			return point_to_segment_distance_batch(p, boundary->p, boundary->num_vertices,ctx->geography);
		}
		return DBL_MAX;
	}else{
		//checking convex
		if(ctx->is_within_query()&&convex_hull){
			double min_dist = DBL_MAX;
			for(int i=0;i<convex_hull->num_vertices-1;i++){
				double dist = point_to_segment_distance(p, convex_hull->p[i], convex_hull->p[i+1], ctx->geography);
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
				return point_to_segment_distance_batch(p, boundary->p, boundary->num_vertices,ctx->geography);
			}
		}else{
			return DBL_MAX;
		}
	}
}


double MyPolygon::distance(MyPolygon *target, query_context *ctx){

	if(raster){

		// both polygons are rasterized and the pixel of the target is bigger
		// then use the target polygon as the host one
		if(target->raster && target->get_rastor()->get_step(false) < get_rastor()->get_step(false)){
			return target->distance(this, ctx);
		}

		double mindist = getMBB()->max_distance(*target->getMBB(), ctx->geography);
		double mbrdist = getMBB()->distance(*target->getMBB(),ctx->geography);
		double min_mbrdist = mbrdist;
		int step = 0;
		double step_size = raster->get_step(ctx->geography);

		vector<Pixel *> needprocess = raster->get_closest_pixels(target->getMBB());
		assert(needprocess.size()>0);
		unsigned short lowx = needprocess[0]->id[0];
		unsigned short highx = needprocess[0]->id[0];
		unsigned short lowy = needprocess[0]->id[1];
		unsigned short highy = needprocess[0]->id[1];
		for(Pixel *p:needprocess){
			lowx = min(lowx, p->id[0]);
			highx = max(highx, p->id[0]);
			lowy = min(lowy, p->id[1]);
			highy = max(highy, p->id[1]);
		}

		while(true){
			struct timeval start = get_cur_time();

			// first of all, expand the circle to involve more pixels
			if(step>0){
				needprocess = raster->expand_radius(lowx,highx,lowy,highy,step);
			}
			ctx->pixel_evaluated.counter += needprocess.size();
			ctx->pixel_evaluated.execution_time += get_time_elapsed(start);

			// all the boxes are scanned (should never happen)
			if(needprocess.size()==0){
				return mindist;
			}

			for(Pixel *cur:needprocess){
				//printf("checking pixel %d %d %d\n",cur->id[0],cur->id[1],cur->status);
				start = get_cur_time();
				ctx->border_evaluated.execution_time += get_time_elapsed(start);
				// note that there is no need to check the edges of
				// this pixel if it is too far from the target
				if(cur->is_boundary() && cur->distance(*target->getMBB(),ctx->geography) < mindist){
					ctx->border_evaluated.counter++;
					// the vector model need be checked.
					start = get_cur_time();
					ctx->border_checked.counter++;
					int edge_num = 0;
					for(edge_range &er:cur->edge_ranges){
						for(int c=er.vstart;c<=er.vend;c++){
							double dist = target->distance(boundary->p[c], ctx);
							mindist = min(dist, mindist);
						}
						edge_num += er.size();
					}
					ctx->edge_checked.counter += cur->num_edges_covered();
					ctx->edge_checked.execution_time += get_time_elapsed(start);
				}
			}
			//log("step:%d #pixels:%ld radius:%f mindist:%f",step, needprocess.size(), mbrdist+step*step_size, mindist);
			needprocess.clear();
			// the minimum distance for now is good enough
			if(mindist<=mbrdist+step*step_size){
				return mindist;
			}
			step++;
		}

		// iterate until the closest pair of edges are found
		assert(false && "happens when there is no boundary pixel, check out the input");
		return DBL_MAX;
	}

	// brute-forcely calculate the distance
	return segment_to_segment_distance_batch(boundary->p, target->boundary->p,
							boundary->num_vertices, target->boundary->num_vertices,ctx->geography);
	return 0;
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
		vector<Pixel *> covered = raster->get_intersect_pixels(b);
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




