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
	size_t max_objects_number = 10000;
	int min_num_vertex = INT_MAX;

	po::options_description desc("query usage");
	desc.add_options()
		("help,h", "produce help message")
		("source,s", po::value<string>(&data_path)->required(), "path to the source")
		("number,n", po::value<size_t>(&max_objects_number), "max number of objects")
		("min_num_vertex,n", po::value<int>(&min_num_vertex), "print the first one with at least this number of vertexes")
		("is_points,p", "the input data are points")
		;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
		cout << desc << "\n";
		exit(0);
	}
	po::notify(vm);
	struct timeval start = get_cur_time();
	if(vm.count("is_points")){
		Point *points;
		size_t point_num = ::load_points_from_path(data_path.c_str(), &points);
		MyMultiPoint mp;
		for(size_t i=0;i<point_num;i++){
			mp.insert(points + i);
		}
		mp.print(max_objects_number);
		logt("%ld points are printed",start, point_num);
	}else{
		query_context ctx;
		size_t num_objects = number_of_objects(data_path.c_str());
		ctx.sample_rate = 1.0*max_objects_number/num_objects;
		vector<MyPolygon *> polygons = load_binary_file(data_path.c_str(), ctx);
		for(MyPolygon *p:polygons){
			if(vm.count("min_num_vertex") && p->get_num_vertices()>=min_num_vertex){
				p->print(false);
				break;
			}
			p->print(false);
		}
		for(MyPolygon *p:polygons){
			delete p;
		}
		logt("%ld polygons are loaded",start, polygons.size());

		polygons.clear();
	}
	return 0;
}
