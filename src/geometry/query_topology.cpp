/*
 * query.cpp
 *
 *  Created on: Jun 2, 2020
 *      Author: teng
 */


#include <float.h>
#include <math.h>
#include <utility>

#include "../include/geometry_computation.h"
#include "../include/MyPolygon.h"

/*
 *
 * containment test
 *
 * */

bool MyPolygon::contain(Point &p){
    return boundary->contain(p);
}

bool contain_rtree(RTNode *node, Point &p, query_context *ctx){
	if(!node->contain(p)){
		return false;
	}
	if(node->is_leaf()){
		assert(node->node_element);
		Point *triangle = (Point *)node->node_element;
		struct timeval start = get_cur_time();
		bool ret = false;
		for(int i=0;i<=2;i++){
			Point *start = triangle+i;
			Point *end = triangle+(i+1)%3;
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
	// todo adjust the lower bound of pixel number when the raster model is usable
	if(raster && get_num_pixels()>5){
		start = get_cur_time();
		int target = raster->get_pixel_id(p);
		auto pix = raster->get_pixels();
		box bx = raster->get_pixel_box(raster->get_x(target), raster->get_y(target));
		double bx_high = bx.high[0];
		if(profile){
			ctx->pixel_evaluated.counter++;
			ctx->pixel_evaluated.execution_time += get_time_elapsed(start);
		}
		if(pix->show_status(target) == IN) {
			return true;
		}
		if(pix->show_status(target) == OUT){
			return false;
		}

		start = get_cur_time();
		bool ret = false;

		// checking the intersection edges in the target pixel
		uint edge_count = 0;
        for(uint16_t e = 0; e < pix->get_num_sequences(target); e ++){    
            auto edges = pix->get_edge_sequence(pix->get_pointer(target) + e);
            auto pos = edges.first;
            for(int k = 0; k < edges.second; k ++){
                int i = pos + k;
                int j = i + 1;  //ATTENTION
                if(((boundary->p[i].y >= p.y) != (boundary->p[j].y >= p.y))){
					double int_x = (boundary->p[j].x - boundary->p[i].x) * (p.y - boundary->p[i].y) / (boundary->p[j].y - boundary->p[i].y) + boundary->p[i].x;
					if(p.x <= int_x && int_x <= bx_high){
						ret = !ret;
					}
				}
            }
			edge_count += edges.second;
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
	}else if(qtree && this->get_num_pixels()>5){
		start = get_cur_time();
		QTNode *tnode = get_qtree()->retrieve(p);
		assert(tnode->isleaf()&&tnode->mbr.contain(p));
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

		// todo for test only, remove in released version
		if(ctx->perform_refine)
		{
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

	if(ctx->use_geos){
		assert(geos_geom && target->geos_geom);
		return geos_geom->covers(target->geos_geom.get());
	}

	if(raster){
		vector<int> pxs = raster->retrieve_pixels(target->getMBB());
		auto pixs = raster->get_pixels();
		auto tpixs = target->raster->get_pixels();
		int etn = 0;
		int itn = 0;
		for(auto p : pxs){
			if(pixs->show_status(p) == OUT){
				etn++;
			}else if(pixs->show_status(p) == IN){
				itn++;
			}
		}
		//log("%d %d %d",etn,itn,pxs.size());
		if(etn == pxs.size()){
			return false;
		}
		if(itn == pxs.size()){
			return true;
		}
		ctx->border_checked.counter++;

		start = get_cur_time();
		if(target->raster){
			vector<pair<int, int>> candidates;
			vector<int> tpxs;
			start = get_cur_time();
			for(auto p : pxs){
				box bx =  raster->get_pixel_box(raster->get_x(p), raster->get_y(p));
				tpxs = target->raster->retrieve_pixels(&bx);
				for(auto p2 : tpxs){
					ctx->pixel_evaluated.counter++;
					// an external pixel of the container intersects an internal
					// pixel of the containee, which means the containment must be false
					if(pixs->show_status(p) == IN) continue;
					if(pixs->show_status(p) == OUT && tpixs->show_status(p2) == IN){
						ctx->pixel_evaluated.execution_time += get_time_elapsed(start,true);
						return false;
					}
					if (pixs->show_status(p) == OUT && tpixs->show_status(p2) == BORDER){
						Point pix_border[5];
						pix_border[0].x = bx.low[0]; pix_border[0].y = bx.low[1];
						pix_border[1].x = bx.low[0]; pix_border[1].y = bx.high[1];
						pix_border[2].x = bx.high[0]; pix_border[2].y = bx.high[1];
						pix_border[3].x = bx.high[0]; pix_border[3].y = bx.low[1];
						pix_border[4].x = bx.low[0]; pix_border[4].y = bx.low[1];
						for (int e = 0; e < tpixs->get_num_sequences(p2); e++){
							auto edges = tpixs->get_edge_sequence(tpixs->get_pointer(p2) + e);
							auto pos = edges.first;
							auto size = edges.second;
							if (segment_intersect_batch(target->boundary->p + pos, pix_border, size, 4, ctx->edge_checked.counter)){
								ctx->edge_checked.execution_time += get_time_elapsed(start, true);
								return false;
							}
						}
					}
					// evaluate the state
					if(pixs->show_status(p) == BORDER && tpixs->show_status(p2) == BORDER){
						candidates.push_back(make_pair(p, p2));
					}
				}
				tpxs.clear();
			}
			ctx->pixel_evaluated.execution_time += get_time_elapsed(start,true);

			for(auto pa : candidates){
				auto p = pa.first;
				auto p2 = pa.second;
				ctx->border_evaluated.counter++;
				// for(edge_range &r:p->edge_ranges){
				// 	for(edge_range &r2:p2->edge_ranges){
				// 		if(segment_intersect_batch(boundary->p+r.vstart, target->boundary->p+r2.vstart, r.size(), r2.size(), ctx->edge_checked.counter)){
				// 			ctx->edge_checked.execution_time += get_time_elapsed(start,true);
				// 			return false;
				// 		}
				// 	}
				// }
				for(int i = 0; i < pixs->get_num_sequences(p); i ++){
					auto r = pixs->get_edge_sequence(pixs->get_pointer(p) + i);
					for(int j = 0; j < tpixs->get_num_sequences(p2); j ++){
						auto r2 = tpixs->get_edge_sequence(tpixs->get_pointer(p2) + j);
						if(segment_intersect_batch(boundary->p+r.first, target->boundary->p+r2.first, r.second, r2.second, ctx->edge_checked.counter)){
							ctx->edge_checked.execution_time += get_time_elapsed(start,true);
							return false;
						}
					}
				}				
			}

			ctx->edge_checked.execution_time += get_time_elapsed(start,true);
		}else{
			vector<int> pxs = raster->retrieve_pixels(target->getMBB());
			vector<int> bpxs;
			for(auto p : pxs){
				if(pixs->show_status(p) == BORDER){
					bpxs.push_back(p);
				}
			}
			for(auto p : pxs){
				for(int i = 0; i < pixs->get_num_sequences(p); i ++){
					auto r = pixs->get_edge_sequence(pixs->get_pointer(p) + i);
					if(segment_intersect_batch(this->boundary->p+r.first, target->boundary->p, r.second, target->boundary->num_vertices, ctx->edge_checked.counter)){
						//logt("%ld boundary %d(%ld) %d(%ld)",start,bpxs.size(),getid(),this->get_num_vertices(),target->getid(), target->get_num_vertices());
						return false;
					}
				}
			}
			return true;
		}
		pxs.clear();
	} else if(qtree) {
		// filtering with the mbr of the target against the qtree
		bool isin = false;
		if(get_qtree()->determine_contain(*(target->getMBB()), isin)){
			return isin;
		}
		ctx->border_checked.counter++;
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
		ctx->border_checked.counter++;

		// use the internal rtree if it is created
		if(rtree){
			for(int i=0;i<target->get_num_vertices();i++){
				if(!contain_rtree(rtree, *target->get_point(i), ctx)){
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
 *
 * intersect
 *
 * */

bool MyPolygon::intersect(MyPolygon *target, query_context *ctx){

	if(!getMBB()->intersect(*target->getMBB())){
		return false;
	}
	ctx->object_checked.counter++;
	if(ctx->use_geos){
		assert(geos_geom && target->geos_geom);
		return geos_geom->intersects(target->geos_geom.get());
	}

	struct timeval start = get_cur_time();

	if(raster){
		vector<int> pxs = raster->retrieve_pixels(target->getMBB());
		auto pixs = raster->get_pixels();
		auto tpixs = target->raster->get_pixels();
		int etn = 0;
		int itn = 0;
		vector<int> bpxs;
		for(auto p : pxs){
			if(pixs->show_status(p) == OUT){
				etn++;
			}else if(pixs->show_status(p) == IN){
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
		ctx->border_checked.counter++;

		start = get_cur_time();
		assert(target->raster);
		vector<pair<int, int>> candidates;
		vector<int> bpxs2;
		start = get_cur_time();
		for(auto p : bpxs){
			box bx = raster->get_pixel_box(raster->get_x(p), raster->get_y(p));
			bpxs2 = target->raster->retrieve_pixels(&bx);
			// cannot determine anything; e-i e-e e-b
			if(pixs->show_status(p) == OUT){
				continue;
			}
			for(auto p2 : bpxs2){
				ctx->pixel_evaluated.counter++;

				// nothing specific; i-e b-e
				if(tpixs->show_status(p2) == OUT){
					continue;
				}

				// must intersect; i-i
				if(pixs->show_status(p) == IN && tpixs->show_status(p2) == OUT){
					ctx->pixel_evaluated.execution_time += get_time_elapsed(start,true);
					return true;
				}

				// b-b b-i i-b
				candidates.push_back(pair<int, int>(p, p2));
			}
			bpxs2.clear();
		}
		ctx->pixel_evaluated.execution_time += get_time_elapsed(start,true);


		for(auto pa : candidates){
			int p = pa.first;
			int p2 = pa.second;
			ctx->border_evaluated.counter++;
			if(pixs->show_status(p) == BORDER && tpixs->show_status(p2) == BORDER){
				for(int i = 0; i < pixs->get_num_sequences(p); i ++){
					auto r = pixs->get_edge_sequence(pixs->get_pointer(p) + i);
					for(int j = 0; j < tpixs->get_num_sequences(p2); j ++){
						auto r2 = tpixs->get_edge_sequence(tpixs->get_pointer(p2) + j);
						if(segment_intersect_batch(boundary->p+r.first, target->boundary->p+r2.first, r.second, r2.second, ctx->edge_checked.counter)){
							ctx->edge_checked.execution_time += get_time_elapsed(start,true);
							return false;
						}
					}
				}				
			} else if (pixs->show_status(p) == BORDER){
				assert(tpixs->show_status(p2) == IN);
				for(int i = 0; i < pixs->get_num_sequences(p); i ++){
					auto r = pixs->get_edge_sequence(pixs->get_pointer(p) + i);
					for(int j = 0; j < r.second;j ++){
                        box bx = target->get_rastor()->get_pixel_box(
                            target->get_rastor()->get_x(p2),
                            target->get_rastor()->get_y(p2));
                        if(bx.intersect(*(boundary->p+r.first+i), *(boundary->p+r.first+i+1))){
							ctx->edge_checked.execution_time += get_time_elapsed(start,true);
							return true;
						}
					}
				}
			} else {
				assert(pixs->show_status(p) == IN);
				for(int i = 0; i < tpixs->get_num_sequences(p2); i ++){
					auto r = tpixs->get_edge_sequence(tpixs->get_pointer(p2) + i);
					for(int j = 0; j < r.second;j ++){
                        box bx = raster->get_pixel_box(raster->get_x(p),
                                                       raster->get_y(p));
                        if(bx.intersect(*(target->boundary->p+r.first+i), *(target->boundary->p+r.first+i+1))){
							ctx->edge_checked.execution_time += get_time_elapsed(start,true);
							return true;
						}
					}
				}
			}
		}

		// when reach here, check if contain
		if(getMBB()->contain(*target->getMBB())){
			Point p(target->getx(0),target->gety(0));
			bool contained = contain(p, ctx,false);
			ctx->edge_checked.execution_time += get_time_elapsed(start,true);
			return contained;
		}

		if(target->getMBB()->contain(*getMBB())){
			Point p(getx(0), gety(0));
			bool contained = target->contain(p, ctx,false);
			ctx->edge_checked.execution_time += get_time_elapsed(start,true);
			return contained;
		}

		ctx->edge_checked.execution_time += get_time_elapsed(start,true);
		// do not intersect
		return false;
	}

	return segment_intersect_batch(this->boundary->p, target->boundary->p, this->boundary->num_vertices, target->boundary->num_vertices, ctx->edge_checked.counter);
}

bool MyPolygon::intersect_box(box *target){
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

/*
 * geos library wrappers
 * */

bool MyPolygon::contain(geos::geom::Geometry *g){
	assert(geos_geom);
	return geos_geom->contains(g);
}
bool MyPolygon::intersect(geos::geom::Geometry *g){
	assert(geos_geom);
	return geos_geom->intersects(g);
}

double MyPolygon::distance(geos::geom::Geometry *g){
	assert(geos_geom);
	return geos_geom->distance(g);
}