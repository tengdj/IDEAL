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
		Point *triangle = (Point *)node->node_element;
		struct timeval start = get_cur_time();
		bool ret = false;
		for(int i=0,j=2;i<=2;j=i++){
			Point *start = triangle+i;
			Point *end = triangle+j;
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

bool MyPolygon::contain(Point &p, query_context *ctx, bool profile){

	// the MBB may not be checked for within query
	if(!mbr->contain(p)){
		return false;
	}
	if(profile){
		ctx->object_checked.counter++;
	}

	struct timeval start = get_cur_time();
	if(raster){
		start = get_cur_time();
		Pixel *target = raster->get_pixel(p);
		if(profile){
			ctx->pixel_evaluated.counter++;
			ctx->pixel_evaluated.execution_time += get_time_elapsed(start);
		}
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
		if(profile){
			ctx->edge_checked.counter += edge_count;
			ctx->edge_checked.execution_time += get_time_elapsed(start);
		}

		// check the crossing nodes on the right bar
		// swap the state of ret if odd number of intersection
		// nodes encountered at the right side of the border
		struct timeval tstart = get_cur_time();
		int nc = raster->count_intersection_nodes(p);
		if(nc%2==1){
			ret = !ret;
		}
		if(profile){
			ctx->intersection_checked.counter += nc;
			ctx->intersection_checked.execution_time += get_time_elapsed(tstart);

			ctx->border_checked.counter++;
			ctx->border_checked.execution_time += get_time_elapsed(start);
			ctx->refine_count++;
		}

		return ret;
	}else if(qtree){
		start = get_cur_time();
		QTNode *tnode = get_qtree()->retrieve(p);
		assert(tnode->isleaf&&tnode->mbr.contain(p));
		if(profile){
			ctx->pixel_evaluated.counter++;
			ctx->pixel_evaluated.execution_time += get_time_elapsed(start);
		}

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
		if(profile){
			ctx->refine_count++;
			ctx->edge_checked.counter += get_num_vertices();
			ctx->edge_checked.execution_time += get_time_elapsed(start);
			ctx->border_checked.counter++;
			ctx->border_checked.execution_time += get_time_elapsed(start);
		}


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
		if(profile){
			ctx->refine_count++;
			ctx->edge_checked.counter += get_num_vertices();
			ctx->edge_checked.execution_time += get_time_elapsed(start);
			ctx->border_checked.counter++;
			ctx->border_checked.execution_time += get_time_elapsed(start);
		}
		return contained;
	}
}

bool MyPolygon::contain(MyPolygon *target, query_context *ctx){
	if(!getMBB()->contain(*target->getMBB())){
		//log("mbb do not contain");
		return false;
	}
	ctx->object_checked.counter++;
	struct timeval start = get_cur_time();

	if(raster){
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
		ctx->pixel_evaluated.counter += pxs.size();
		ctx->pixel_evaluated.execution_time += get_time_elapsed(start,true);
		if(etn == pxs.size()){
			return false;
		}
		if(itn == pxs.size()){
			return true;
		}
		ctx->border_checked.counter++;

		if(target->raster){
			vector<Pixel *> bpxs2;
			start = get_cur_time();
			for(Pixel *p:bpxs){
				bpxs2 = target->raster->retrieve_pixels(p);
				for(Pixel *p2:bpxs2){
					for(edge_range &r:p->edge_ranges){
						for(edge_range &r2:p2->edge_ranges){
							if(segment_intersect_batch(boundary->p+r.vstart, target->boundary->p+r2.vstart, r.size(), r2.size(), ctx->edge_checked.counter)){
								//ctx->edge_checked.execution_time += get_time_elapsed(start,true);
								return false;
							}
						}
					}
				}
				bpxs2.clear();
			}
			//ctx->edge_checked.execution_time += get_time_elapsed(start,true);
		}else{
			for(Pixel *p:bpxs){
				for(edge_range &r:p->edge_ranges){
					if(segment_intersect_batch(this->boundary->p+r.vstart, target->boundary->p, r.size(), target->boundary->num_vertices, ctx->edge_checked.counter)){
						//logt("%ld boundary %d(%ld) %d(%ld)",start,bpxs.size(),getid(),this->get_num_vertices(),target->getid(), target->get_num_vertices());
						return false;
					}
				}
			}
		}
		bpxs.clear();
	} else if(qtree) {
		// filtering with the mbr of the target against the qtree
		bool isin = false;
		if(get_qtree()->determine_contain(*(target->getMBB()), isin)){
			return isin;
		}
		// otherwise, checking all the edges to make sure no intersection
		if(segment_intersect_batch(boundary->p, target->boundary->p, boundary->num_vertices, target->boundary->num_vertices, ctx->edge_checked.counter)){
			return false;
		}
	} else {
		// filtering with mer, mbr, and convex hull
		if(mer){
			// filter with the mer of source and mbr of the target
			if(mer->contain(*target->getMBB())){
				return true;
			}

			// filter with the convex hull of target
			if(target->convex_hull){
				Point mer_vertices[5];
				mer->to_array(mer_vertices);
				if(!segment_intersect_batch(mer_vertices, target->convex_hull->p, 5, target->convex_hull->num_vertices, ctx->edge_checked.counter)){
					if(mer->contain(convex_hull->p[0])){
						return true;
					}
				}
			}
		}

		// further filtering with the mbr of the target
		Point mbb_vertices[5];
		target->mbr->to_array(mbb_vertices);
		// no intersection between this polygon and the mbr of the target polygon
		if(!segment_intersect_batch(boundary->p, mbb_vertices, boundary->num_vertices, 5, ctx->edge_checked.counter)){
			// the target must be the one which is contained (not contain) as its mbr is contained
			if(contain(mbb_vertices[0], ctx)){
				return true;
			}
		}

		// when reach here, we have no choice but evaluate all edge pairs

		// use rtree if it is created
		if(rtree){
			for(int i=0;i<target->get_num_vertices();i++){
				if(!rtree->contain(*target->get_point(i))){
					return false;
				}
			}
			return true;
		}

		// otherwise, checking all the edges to make sure no intersection
		if(segment_intersect_batch(boundary->p, target->boundary->p, boundary->num_vertices, target->boundary->num_vertices, ctx->edge_checked.counter)){
			return false;
		}
	}

	// this is the last step for all the cases, when no intersection segment is identified
	// pick one point from the target and it must be contained by this polygon
	Point p(target->getx(0),target->gety(0));
	return contain(p, ctx,false);

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
					Point *triangle = (Point *)pr.first->node_element;
					for(int i=0,j=2;i<=2;j=i++){
						double dist = point_to_segment_distance(p, triangle[i], triangle[j],ctx->geography);
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
// calculate the distance with rtree from a segment
double MyPolygon::distance_rtree(Point &start, Point &end, query_context *ctx){
	assert(rtree);

	// start a breadth first search
	queue<RTNode *> pixq;
	for(RTNode *p:rtree->children){
		pixq.push(p);
	}
	// set this value as the MINMAXDIST
	double mindist = DBL_MAX;
	vector<std::pair<RTNode *, double>> candidates;
	while(pixq.size()>0){
		for(int i=0;i<pixq.size();i++){
			RTNode *pix = pixq.front();
			pixq.pop();
			// this node and its descendants must not be candidate
			double dist = pix->distance(start, end, ctx->geography);
			if(dist >= mindist){
				continue;
			}
			// maximum possible distance between the target and this node
			double maxdist = pix->max_distance(start, end, ctx->geography);
			mindist = min(maxdist, mindist);
			candidates.push_back(std::pair<RTNode *, double>(pix, dist));
		}

		for(std::pair<RTNode *, double> pr:candidates){
			// minimum possible distance of this node
			// must be smaller than the current minimum
			if(pr.second<=mindist){
				// is leaf
				if(pr.first->is_leaf()){
					Point *triangle = (Point *)pr.first->node_element;
					for(int i=0;i<=2;i++){
						double dist = segment_to_segment_distance(start, end, triangle[i], triangle[(i+1)%3],ctx->geography);
						mindist = min(mindist, dist);
						// close enough for a within query
						if(ctx->within(mindist)){
							return mindist;
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
		candidates.clear();
	}

	return mindist;
}

double MyPolygon::distance(Point &p, query_context *ctx, bool profile){

	// distance is 0 if contained by the polygon
	double mindist = getMBB()->max_distance(p, ctx->geography);

	struct timeval query_start = get_cur_time();
	bool contained = contain(p, ctx, profile);
	if(profile){
		ctx->contain_check.execution_time += get_time_elapsed(query_start,true);
		ctx->contain_check.counter++;
	}

	if(contained){
		return 0;
	}
	if(profile){
		ctx->object_checked.counter++;
	}

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
				if(profile){
					ctx->refine_count++;
				}
				return point_to_segment_distance_batch(p, boundary->p, boundary->num_vertices,ctx->geography);
			}
			//if(profile)
			{
				ctx->pixel_evaluated.counter += needprocess.size();
				ctx->pixel_evaluated.execution_time += get_time_elapsed(start, true);
			}

			for(Pixel *cur:needprocess){
				//printf("checking pixel %d %d %d\n",cur->id[0],cur->id[1],cur->status);
				if(cur->is_boundary()){
					there_is_border = true;
					start = get_cur_time();
					//if(profile)
					{
						ctx->border_evaluated.counter++;
					}
					// no need to check the edges of this pixel
					double mbr_dist = cur->distance(p, ctx->geography);
					if(profile){
						ctx->border_evaluated.execution_time += get_time_elapsed(start);
					}
					// skip the boundary pixels that is further than the
					// current minimum
					if(mbr_dist >= mindist){
						continue;
					}
					border_checked = true;

					// the vector model need be checked.
					start = get_cur_time();
					if(profile){
						ctx->border_checked.counter++;
					}

					for(edge_range &rg:cur->edge_ranges){
						for (int i = rg.vstart; i <= rg.vend; i++) {
							ctx->edge_checked.counter ++;
							double dist = point_to_segment_distance(p, *get_point(i), *get_point(i+1),ctx->geography);
							mindist = min(mindist, dist);
							if(ctx->within(mindist)){
								ctx->edge_checked.execution_time += get_time_elapsed(start);
								return dist;
							}

						}
					}
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
		if(profile){
			ctx->refine_count++;
		}
		// IDEAL return
		return mindist;

	}else if(qtree){
		if(ctx->is_within_query()&&
				!qtree->within(p, ctx->within_distance)){
			return DBL_MAX;
		}
		if(profile){
			ctx->refine_count++;
		}
		if(ctx->perform_refine){
			if(profile){
				ctx->edge_checked.counter += this->get_num_vertices();
			}
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
			if(min_dist>ctx->within_distance){
				return min_dist;
			}
		}
		if(profile){
			ctx->refine_count++;
		}
		//SIMPVEC return
		if(ctx->perform_refine){
			if(rtree){
				return distance_rtree(p,ctx);
			}else{
				if(profile){
					ctx->edge_checked.counter += get_num_vertices();
				}
				return point_to_segment_distance_batch(p, boundary->p, boundary->num_vertices,ctx->geography);
			}
		}else{
			return DBL_MAX;
		}
	}
}

// get the distance from pixel pix to polygon target
double MyPolygon::distance(MyPolygon *target, Pixel *pix, query_context *ctx, bool profile){
	double mindist = DBL_MAX;
	assert(target->raster);
	assert(pix->is_boundary());

	if(raster){
		mindist = getMBB()->max_distance(*pix, ctx->geography);
		const double mbrdist = getMBB()->distance(*pix,ctx->geography);
		double min_mbrdist = mbrdist;
		int step = 0;
		double step_size = raster->get_step(ctx->geography);
		// initialize the seed closest pixels
		vector<Pixel *> needprocess = raster->get_closest_pixels(pix);
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

			// for later steps, expand the circle to involve more pixels
			if(step>0){
				needprocess = raster->expand_radius(lowx,highx,lowy,highy,step);
			}
			//ctx->pixel_evaluated.counter += needprocess.size();
			//ctx->pixel_evaluated.execution_time += get_time_elapsed(start, true);

			// all the boxes are scanned (should never happen)
			if(needprocess.size()==0){
				return mindist;
			}

			for(Pixel *cur:needprocess){
				if(profile){
					ctx->pixel_evaluated.counter++;
				}
				//printf("checking pixel %d %d %d\n",cur->id[0],cur->id[1],cur->status);
				// note that there is no need to check the edges of
				// this pixel if it is too far from the target
				if(cur->is_boundary()){
					start = get_cur_time();
					bool toofar = (cur->distance(*pix,ctx->geography) >= mindist);
					if(profile){
						ctx->border_evaluated.counter++;
						ctx->border_evaluated.execution_time += get_time_elapsed(start, true);
					}
					if(toofar){
						continue;
					}
					//ctx->border_evaluated.counter++;
					// the vector model need be checked.
					start = get_cur_time();
					if(profile){
						ctx->border_checked.counter++;
					}
					for(edge_range &pix_er:pix->edge_ranges){
						for(edge_range &cur_er:cur->edge_ranges){
							double dist;
							if(ctx->is_within_query()){
								dist = segment_to_segment_within_batch(target->boundary->p+pix_er.vstart,
													boundary->p+cur_er.vstart, pix_er.size(), cur_er.size(),
													ctx->within_distance, ctx->geography, ctx->edge_checked.counter);
							}else{
								dist = segment_to_segment_distance_batch(target->boundary->p+pix_er.vstart,
													boundary->p+cur_er.vstart, pix_er.size(), cur_er.size(), ctx->geography);
								if(profile){
									ctx->edge_checked.counter += pix_er.size()*cur_er.size();
								}
							}
							mindist = min(dist, mindist);
							if(ctx->within(mindist)){
								return mindist;
							}
						}
					}
					if(profile){
						ctx->edge_checked.execution_time += get_time_elapsed(start, true);
					}
				}
			}
			//log("step:%d #pixels:%ld radius:%f mindist:%f",step, needprocess.size(), mbrdist+step*step_size, mindist);
			needprocess.clear();
			double min_possible = mbrdist+step*step_size;
			// the minimum distance for now is good enough for three reasons:
			// 1. current minimum distance is smaller than any further distance
			// 2. for within query, current minimum is close enough
			// 3. for within query, current minimum could never be smaller than the threshold
			if(mindist <= min_possible
			   || (ctx->within(mindist))
			   || (ctx->is_within_query() && ctx->within_distance < min_possible)){
				return mindist;
			}
			step++;
		}

	}else{
		for(edge_range &er:pix->edge_ranges){
			double dist;
			if(ctx->is_within_query()){
				dist = segment_to_segment_within_batch(target->boundary->p+er.vstart,
									boundary->p, er.size(), boundary->num_vertices, ctx->within_distance, ctx->geography, ctx->edge_checked.counter);
			}else{
				dist = segment_to_segment_distance_batch(target->boundary->p+er.vstart,
									boundary->p, er.size(), boundary->num_vertices, ctx->geography);
				if(profile){
					ctx->edge_checked.counter += er.size()*boundary->num_vertices;
				}
			}

			mindist = min(dist, mindist);
			if(ctx->within(mindist)){
				return mindist;
			}
		}
		return mindist;
	}

	return mindist;
}

double MyPolygon::distance(MyPolygon *target, query_context *ctx){

	ctx->object_checked.counter++;

	if(raster){
		// both polygons are rasterized and the pixel of the target is larger
		// then swap the role of source and target, and use the target as the host one
		if(target->raster && target->get_rastor()->get_step(false) < get_rastor()->get_step(false)){
			ctx->object_checked.counter--;
			return target->distance(this, ctx);
		}

		double mindist = getMBB()->max_distance(*target->getMBB(), ctx->geography);
		const double mbrdist = getMBB()->distance(*target->getMBB(),ctx->geography);
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
			//ctx->pixel_evaluated.counter += needprocess.size();
			//ctx->pixel_evaluated.execution_time += get_time_elapsed(start, true);

			// all the boxes are scanned (should never happen)
			if(needprocess.size()==0){
				return mindist;
			}

			for(Pixel *cur:needprocess){
				//printf("checking pixel %d %d %d\n",cur->id[0],cur->id[1],cur->status);
				// note that there is no need to check the edges of
				// this pixel if it is too far from the target
				if(cur->is_boundary() && cur->distance(*target->getMBB(),ctx->geography) < mindist){
					// the vector model need be checked.
					// do a polygon--pixel distance calculation
					double dist = target->distance(this, cur, ctx, true);
					mindist = min(dist, mindist);
					if(ctx->within(mindist)){
						return mindist;
					}
				}
			}
			//log("step:%d #pixels:%ld radius:%f mindist:%f",step, needprocess.size(), mbrdist+step*step_size, mindist);
			needprocess.clear();
			double min_possible = mbrdist+step*step_size;
			// the minimum distance for now is good enough for three reasons:
			// 1. current minimum distance is smaller than any further distance
			// 2. for within query, current minimum is close enough
			// 3. for within query, current minimum could never be smaller than the threshold
			if(mindist <= min_possible
					|| (ctx->within(mindist))
					|| (ctx->is_within_query() && ctx->within_distance < min_possible)){
				return mindist;
			}
			step++;
		}

		// iterate until the closest pair of edges are found
		assert(false && "happens when there is no boundary pixel, check out the input");
		return DBL_MAX;
	}else if(qtree){

		// checking the qtree for filtering
		if(ctx->is_within_query()){
			for(int i=0;i<target->boundary->num_vertices-1;i++){
				// if any edge is within the distance after checking the QTree, get the exact distance from it to the source polygon
				if(qtree->within(target->boundary->p[i], target->boundary->p[i+1], ctx->within_distance)){
					double dist = segment_to_segment_distance_batch(boundary->p, target->boundary->p+i,
							boundary->num_vertices, 1,ctx->geography);
					if(dist <= ctx->within_distance){
						return dist;
					}
				}
			}
		}
		// qtree do not support general distance calculation, do computation when failed filtering
		return segment_to_segment_distance_batch(boundary->p, target->boundary->p,
								boundary->num_vertices, target->boundary->num_vertices,ctx->geography);
	}else{
		//checking convex for filtering
		if(ctx->is_within_query() && convex_hull && target->convex_hull){
			double dist = segment_to_segment_distance_batch(convex_hull->p, target->convex_hull->p,
					convex_hull->num_vertices, target->convex_hull->num_vertices, ctx->geography);
			// the convex_hull is a conservative approximation of the original polygon,
			// so the distance between the convex hulls of two polygon is smaller than
			// the real distance
			if(dist>ctx->within_distance){
				log("%f",dist);
				return dist;
			}
		}

		//SIMPVEC return
		if(rtree){
			double mindist = DBL_MAX;
			for(int i=0;i<target->boundary->num_vertices-1;i++){
				double dist = distance_rtree(target->boundary->p[i], target->boundary->p[i+1], ctx);
				mindist = min(mindist, dist);
				//log("%f",mindist);

				if(ctx->within(mindist)){
					return mindist;
				}
			}
			return mindist;
		}else{
			return segment_to_segment_distance_batch(boundary->p, target->boundary->p,
									boundary->num_vertices, target->boundary->num_vertices,ctx->geography);
		}
	}
	assert(false);
	// should never reach here
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




