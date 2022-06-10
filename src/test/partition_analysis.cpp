/*
 * partition_analysis.cpp
 *
 *  Created on: Jun 8, 2022
 *      Author: teng
 */

#include "partition.h"


bool assign(Tile *tile, void *arg){
	tile->insert((box *)arg, (box *)arg);
	return true;
}

bool search(Tile *tile, void *arg){
	query_context *ctx = (query_context *)arg;
	ctx->found += tile->lookup_count((Point *)ctx->target);
	return true;
}

int main(int argc, char** argv) {

	const char *mbr_path = argv[1];
	const char *point_path = argv[2];
	int partition_num = atoi(argv[3]);
	PARTITION_TYPE ptype = parse_partition_type(argv[4]);

	struct timeval start = get_cur_time();
	box *boxes;
	size_t box_num = load_boxes_from_file(mbr_path, &boxes);
	logt("%ld objects are loaded",start, box_num);

	Point *points;
	size_t points_num = load_points_from_path(point_path, &points);
	logt("%ld points are loaded", start, points_num);
	for(size_t s=64;s<=64;s*=2){
		double sample_rate = 0.1*s/64;

		vector<box *> objects;
		for(size_t i=0;i<box_num;i++){
			if(tryluck(sample_rate)){
				objects.push_back(boxes+i);
			}
		}
		logt("%ld objects are sampled",start, objects.size());

		vector<Point *> targets;
		for(size_t i=0;i<points_num;i++){
			if(tryluck(sample_rate)){
				targets.push_back(points+i);
			}
		}
		logt("%ld targets are sampled",start, targets.size());

		for(size_t pr=1;pr<=16;pr*=2){
			size_t true_partition_num = partition_num*pr;
			vector<Tile *> tiles = genschema(objects, true_partition_num, ptype);
			RTree<Tile *, double, 2, double> tree;

			for(Tile *t:tiles){
				tree.Insert(t->low, t->high, t);
			}
			logt("%ld tiles are generated",start, tiles.size());

			for(box *obj:objects){
				//obj->print();
				tree.Search(obj->low, obj->high, assign, (void *)obj);
			}

			size_t total_num = 0;
			for(Tile *tile:tiles){
				total_num += tile->get_objnum();
			}
			logt("partitioning data",start);
			//print_boxes(parts);
			query_context ctx;
			for(Point *p:targets){
				ctx.target = (void *)p;
				tree.Search((double *)p, (double *)p, search, (void *)&ctx);
			}
			logt("%.4f boundary %.4f average hit", start, 100.0*(total_num-objects.size())/objects.size(), 1.0*ctx.found/targets.size());

			for(Tile *tile:tiles){
				delete tile;
			}
			tiles.clear();
		}
		objects.clear();
		targets.clear();
	}


	delete []boxes;
	delete []points;
	return 0;
}
