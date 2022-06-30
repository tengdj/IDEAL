/*
 * subset.cpp
 *
 *  Created on: Jun 8, 2022
 *      Author: teng
 */

#include "partition.h"
#include <boost/sort/sort.hpp>

namespace po = boost::program_options;

bool compare(MyPolygon *a, MyPolygon *b){
	return a->getMBB()->area()>b->getMBB()->area();
}

int main(int argc, char** argv) {

	string in_path;
	string out_path;


	po::options_description desc("query usage");
	desc.add_options()
		("help,h", "produce help message")
		("input,i", po::value<string>(&in_path)->required(), "path to the source")
		("output,o", po::value<string>(&out_path)->required(), "path to the target")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
		cout << desc << "\n";
		exit(0);
	}
	po::notify(vm);

	query_context ctx;
	vector<MyPolygon *> polygons = load_binary_file(in_path.c_str(), ctx);
	struct timeval start = get_cur_time();

	boost::sort::block_indirect_sort(polygons.begin(), polygons.end(), compare);
	logt("sorting",start, polygons.size());
#pragma omp parallel for
	for(size_t i=0;i<polygons.size()/4;i++){
		delete polygons[i];
	}
	polygons.erase(polygons.begin(), polygons.begin()+polygons.size()/4);
	logt("selected %ld polygons",start, polygons.size());

	dump_polygons_to_file(polygons, out_path.c_str());
	logt("stored %ld polygons to %s",start, polygons.size(), out_path.c_str());
#pragma omp parallel for
	for(MyPolygon *p:polygons){
		delete p;
	}
	polygons.clear();

	return 0;
}
