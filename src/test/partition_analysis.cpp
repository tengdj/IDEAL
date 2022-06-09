/*
 * partition_analysis.cpp
 *
 *  Created on: Jun 8, 2022
 *      Author: teng
 */

#include "partition.h"




int main(int argc, char** argv) {
	query_context global_ctx;

	box *boxes;
	size_t box_num = load_boxes_from_file(argv[1], &boxes);
	vector<box *> mbrs;
	for(size_t i=0;i<box_num;i++){
		mbrs.push_back(boxes+i);
	}
	//double *points;
	//size_t points_num = load_points_from_path(argv[2], &points);

	vector<box *> parts = genschema(mbrs, atoi(argv[3]), parse_partition_type(argv[4]));

	print_boxes(parts);
	mbrs.clear();
	delete []boxes;
	return 0;
}
