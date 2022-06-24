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
box china(73.47431,3.393568,134.86609,53.65586);

int main(int argc, char** argv) {

	string data_path;
	size_t max_points_number = 10000;

	po::options_description desc("query usage");
	desc.add_options()
		("help,h", "produce help message")
		("source,s", po::value<string>(&data_path)->required(), "path to the source")
		("points_number,n", po::value<size_t>(&max_points_number), "max number of points")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
		cout << desc << "\n";
		exit(0);
	}
	po::notify(vm);
	struct timeval start = get_cur_time();
	Point *points;
	size_t point_num = ::load_points_from_path(data_path.c_str(), &points);
	MyMultiPoint mp;
	for(size_t i=0;i<point_num;i++){
		mp.insert(points + i);
	}
	logt("%ld points are loaded",start, point_num);
	mp.print(max_points_number);
	return 0;
}
