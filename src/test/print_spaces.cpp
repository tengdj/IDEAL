/*
 * partition_analysis.cpp
 *
 *  Created on: Jun 8, 2022
 *      Author: teng
 */

#include "MyPolygon.h"
#include "query_context.h"

namespace po = boost::program_options;

// functions for sampling
int main(int argc, char** argv) {

	vector<string> files;
	::list_files(argv[1], files);
	vector<box *> boxes;
	for(string &s:files){
		box b = ::universe_space(s.c_str());
		boxes.push_back(new box(b));
	}
	print_boxes(boxes);

	for(box *b:boxes){
		delete b;
	}
	boxes.clear();
	return 0;
}
