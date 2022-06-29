/*
 * partition_analysis.cpp
 *
 *  Created on: Jun 8, 2022
 *      Author: teng
 */

#include "partition.h"

namespace po = boost::program_options;

int main(int argc, char** argv) {

	string in_path;
	string out_path;

	double sample_rate = 0.01;

	po::options_description desc("query usage");
	desc.add_options()
		("help,h", "produce help message")
		("is_point,p", "the input and output are points")
		("input,i", po::value<string>(&in_path)->required(), "path to the source")
		("output,o", po::value<string>(&out_path)->required(), "path to the target")

		("sample_rate,r", po::value<double>(&sample_rate), "the sample rate (0.01 by default)")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
		cout << desc << "\n";
		exit(0);
	}
	po::notify(vm);

	struct timeval start = get_cur_time();

	if(vm.count("is_point")){
		Point *points;
		size_t points_num = load_points_from_path(in_path.c_str(), &points);
		vector<Point> sampled;
		for(size_t i=0;i<points_num;i++){
			if(tryluck(sample_rate)){
				sampled.push_back(points[i]);
			}
			if(i%100==0){
				log_refresh("sampled %.2f\%",100.0*i/points_num);
			}
		}

		dump_to_file(out_path.c_str(), (char *)&sampled[0], sizeof(Point)*sampled.size());
		log("%ld points are sampled", sampled.size());
		sampled.clear();
		delete []points;
	}else{
		query_context ctx;
		ctx.sample_rate = sample_rate;
		vector<MyPolygon *> polygons = load_binary_file(in_path.c_str(), ctx);
		dump_polygons_to_file(polygons, out_path.c_str());
		for(MyPolygon *p:polygons){
			delete p;
		}
		polygons.clear();
	}
	return 0;
}
