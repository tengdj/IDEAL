/*
 * partition_analysis.cpp
 *
 *  Created on: Jun 8, 2022
 *      Author: teng
 */

#include "partition.h"




int main(int argc, char** argv) {
	query_context global_ctx;
	vector<Pixel *> mbrs;

	vector<MyPolygon *> polygons = MyPolygon::load_binary_file(argv[1], global_ctx);
	for(MyPolygon *p:polygons){
		mbrs.push_back(p->getMBB());
	}
	//double *points;
	//size_t points_num = load_points_from_path(argv[2], &points);


	vector<Pixel *> parts = genschema_qt(mbrs, atoi(argv[3]));
	print_boxes(parts);

	return 0;
}
