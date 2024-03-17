/*
 * query_distance.cpp
 *
 *  Created on: Apr 25, 2023
 *      Author: teng
 */

#include "../include/MyPolygon.h"

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
					for(int i=0;i<=2;i++){
						double dist = point_to_segment_distance(p, triangle[i], triangle[(i+1)%3],ctx->geography);
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

#ifdef USE_GPU
	if(ctx->gpu){
		return point_to_segment_sequence_distance_gpu(p, boundary->p, boundary->num_vertices,ctx->geography);
	}
#endif

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
		int closest = raster->get_closest_pixel(p);
		int step = 0;
		double step_size = raster->get_step(ctx->geography);
		vector<int> needprocess;

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
				assert(false&&"should not evaluated all boxes");
				if(profile){
					ctx->refine_count++;
				}
				return boundary->distance(p, ctx->geography);
			}
			if(profile)	{
				ctx->pixel_evaluated.counter += needprocess.size();
				ctx->pixel_evaluated.execution_time += get_time_elapsed(start, true);
			}

			auto pixs = raster->get_pixels();
			for(auto cur : needprocess){
				//printf("checking pixel %d %d %d\n",cur->id[0],cur->id[1],cur->status);
				if(pixs->show_status(cur) == BORDER){
					start = get_cur_time();
					if(profile)
					{
						ctx->border_evaluated.counter++;
					}
					// no need to check the edges of this pixel
					if(profile){
						ctx->border_evaluated.execution_time += get_time_elapsed(start);
					}

					box cur_box = raster->get_pixel_box(raster->get_x(cur), raster->get_y(cur));
					double mbr_dist = cur_box.distance(p, ctx->geography);
					// skip the pixels that is further than the current minimum
					if(mbr_dist >= mindist){
						continue;
					}

					// the vector model need be checked.
					start = get_cur_time();
					if(profile){
						ctx->border_checked.counter++;
					}

					for(int i = 0; i < raster->get_num_sequences(cur); i ++){
						auto rg = pixs->get_edge_sequence(pixs->get_pointer(cur) + i);
						for(int j = 0; j < rg.second; j ++){
							auto r = rg.first + j;
							ctx->edge_checked.counter ++;
							double dist = point_to_segment_distance(p, *get_point(r), *get_point(r+1), ctx->geography);
							mindist = min(mindist, dist);
						}
					}
					ctx->edge_checked.execution_time += get_time_elapsed(start);
				}
			}
			//printf("point to polygon distance - step:%d #pixels:%ld radius:%f mindist:%f\n",step,needprocess.size(),mbrdist+step*step_size,mindist);
			needprocess.clear();

			// for within query, return if the current minimum is close enough
			if(ctx->within(mindist)){
				return mindist;
			}
			step++;
			double minrasterdist = raster->get_possible_min(p, closest, step, ctx->geography);
			//cout<<step<<" "<<mindist<<" "<<minrasterdist<<endl;

			// close enough
			if(mindist < minrasterdist){
				break;
			}
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
			return boundary->distance(p, ctx->geography);
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
				return boundary->distance(p, ctx->geography);
			}
		}else{
			return DBL_MAX;
		}
	}
}

// get the distance from pixel pix to polygon target
double MyPolygon::distance(MyPolygon *target, int pix, query_context *ctx, bool profile){
	assert(target->raster);
	assert(raster);
	auto r_pixs = raster->get_pixels();
    auto t_pixs = target->raster->get_pixels();
	assert(r_pixs->show_status(pix) == BORDER);
	struct timeval start = get_cur_time();

	auto pix_x = raster->get_x(pix);
    auto pix_y = raster->get_y(pix);
    auto pix_box = raster->get_pixel_box(pix_x, pix_y);
	double mindist = getMBB()->max_distance(pix_box, ctx->geography);
	double mbrdist = getMBB()->distance(pix_box, ctx->geography);
	double min_mbrdist = mbrdist;
	int step = 0;
	double step_size = raster->get_step(ctx->geography);
	// initialize the seed closest pixels
	vector<int> needprocess = target->raster->get_closest_pixels(pix_box);
	assert(needprocess.size()>0);
	unsigned short lowx = target->raster->get_x(needprocess[0]);
	unsigned short highx = target->raster->get_x(needprocess[0]);
	unsigned short lowy = target->raster->get_y(needprocess[0]);
	unsigned short highy = target->raster->get_y(needprocess[0]);
	for(auto p : needprocess){
		lowx = min(lowx, (unsigned short)target->raster->get_x(p));
		highx = max(highx, (unsigned short)target->raster->get_x(p));
		lowy = min(lowy, (unsigned short)target->raster->get_y(p));
		highy = max(highy, (unsigned short)target->raster->get_y(p));
	}

	ctx->pixel_evaluated.execution_time += get_time_elapsed(start);

	while(true){
		struct timeval start = get_cur_time();

		// for later steps, expand the circle to involve more pixels
		if(step>0){
			needprocess = target->raster->expand_radius(lowx,highx,lowy,highy,step);
		}
		ctx->pixel_evaluated.counter += needprocess.size();
		ctx->pixel_evaluated.execution_time += get_time_elapsed(start, true);

		// all the boxes are scanned (should never happen)
		if(needprocess.size()==0){
			return mindist;
		}

		for(auto cur : needprocess){
			if(profile){
				ctx->pixel_evaluated.counter++;
			}
			//printf("checking pixel %d %d %d\n",cur->id[0],cur->id[1],cur->status);
			// note that there is no need to check the edges of
			// this pixel if it is too far from the target
			auto cur_x = target->raster->get_x(cur);
			auto cur_y = target->raster->get_y(cur);

			if(t_pixs->show_status(cur) == BORDER){
				start = get_cur_time();
				bool toofar = (target->raster->get_pixel_box(cur_x, cur_y).distance(pix_box,ctx->geography) >= mindist);
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
				for(int i = 0; i < raster->get_num_sequences(pix); i ++){
					auto pix_er = r_pixs->get_edge_sequence(r_pixs->get_pointer(pix) + i);
					for(int j = 0; j < target->raster->get_num_sequences(cur); j ++){
						auto cur_er = t_pixs->get_edge_sequence(t_pixs->get_pointer(cur) + j);
						if(cur_er.second < 2 || pix_er.second < 2) continue;
						double dist;
						if(ctx->is_within_query()){
							dist = segment_to_segment_within_batch(target->boundary->p+cur_er.first,
												boundary->p+pix_er.first, cur_er.second, pix_er.second,
												ctx->within_distance, ctx->geography, ctx->edge_checked.counter);
						}else{
							dist = segment_sequence_distance(target->boundary->p+cur_er.first,
												boundary->p+pix_er.first, cur_er.second, pix_er.second, ctx->geography);
							if(profile){
								ctx->edge_checked.counter += pix_er.second*cur_er.second;
							}
						}
						mindist = min(dist, mindist);
						if(ctx->within(mindist)){
							if(profile){
								ctx->edge_checked.execution_time += get_time_elapsed(start, true);
							}
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

	return mindist;
}

double MyPolygon::distance(MyPolygon *target, query_context *ctx){

	ctx->object_checked.counter++;

	if(ctx->use_geos){
		assert(geos_geom && target->geos_geom);
		return geos_geom->distance(target->geos_geom.get());
	}

	if(raster){
		timeval start = get_cur_time();
		// both polygons are rasterized and the pixel of the target is larger
		// then swap the role of source and target, and use the target as the host one
		// currently we put this optimization in the caller
		if(target->raster && target->get_rastor()->get_step(false) > get_rastor()->get_step(false)){
			ctx->object_checked.counter--;
			ctx->pixel_evaluated.execution_time += get_time_elapsed(start);
			return target->distance(this, ctx);
		}

		double mindist = getMBB()->max_distance(*target->getMBB(), ctx->geography);
		const double mbrdist = getMBB()->distance(*target->getMBB(),ctx->geography);
		double min_mbrdist = mbrdist;
		int step = 0;
		double step_size = raster->get_step(ctx->geography);
		auto r_pixs = raster->get_pixels();
		auto t_pixs = target->raster->get_pixels();

		vector<int> needprocess = raster->get_closest_pixels(*target->getMBB());
		assert(needprocess.size()>0);
		unsigned short lowx = raster->get_x(needprocess[0]);
		unsigned short highx = raster->get_x(needprocess[0]);
		unsigned short lowy = raster->get_y(needprocess[0]);
		unsigned short highy = raster->get_y(needprocess[0]);
		for(auto p : needprocess){
			lowx = min(lowx, (unsigned short)raster->get_x(p));
			highx = max(highx, (unsigned short)raster->get_x(p));
			lowy = min(lowy, (unsigned short)raster->get_y(p));
			highy = max(highy, (unsigned short)raster->get_y(p));
		}
		ctx->pixel_evaluated.execution_time += get_time_elapsed(start);

		while(true){
			struct timeval start = get_cur_time();

			// first of all, expand the circle to involve more pixels
			if(step>0){
				needprocess = raster->expand_radius(lowx,highx,lowy,highy,step);
			}
			ctx->pixel_evaluated.counter += needprocess.size();
			ctx->pixel_evaluated.execution_time += get_time_elapsed(start, true);

			// all the boxes are scanned (should never happen)
			if(needprocess.size()==0){
				return mindist;
			}

			for(auto cur : needprocess){
				//printf("checking pixel %d %d %d\n",cur->id[0],cur->id[1],cur->status);
				// note that there is no need to check the edges of
				// this pixel if it is too far from the target
				auto cur_x = raster->get_x(cur);
				auto cur_y = raster->get_y(cur); 
				if(r_pixs->show_status(cur) == BORDER && raster->get_pixel_box(cur_x, cur_y).distance(*target->getMBB(),ctx->geography) < mindist){
					// the vector model need be checked.
					// do a polygon--pixel distance calculation
					double dist = distance(target, cur, ctx, true);
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
					double dist = segment_sequence_distance(boundary->p, target->boundary->p+i,
							boundary->num_vertices, 1,ctx->geography);
					if(dist <= ctx->within_distance){
						return dist;
					}
				}
			}
		}

		// qtree do not support general distance calculation, do computation when failed filtering
		return segment_sequence_distance(boundary->p, target->boundary->p,
								boundary->num_vertices, target->boundary->num_vertices,ctx->geography);
	}else{
		//checking convex for filtering
		if(ctx->is_within_query() && convex_hull && target->convex_hull){
			double dist = segment_sequence_distance(convex_hull->p, target->convex_hull->p,
					convex_hull->num_vertices, target->convex_hull->num_vertices, ctx->geography);
			// the convex_hull is a conservative approximation of the original polygon,
			// so the distance between the convex hulls of two polygon is smaller than
			// the real distance
			if(dist>ctx->within_distance){
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
			return segment_sequence_distance(boundary->p, target->boundary->p,
									boundary->num_vertices, target->boundary->num_vertices,ctx->geography);
		}
	}
	assert(false);
	// should never reach here
}
