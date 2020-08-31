/*
 * Parser.cpp
 *
 *  Created on: May 9, 2020
 *      Author: teng
 */

#include "../geometry/MyPolygon.h"
#include <fstream>
#include <boost/program_options.hpp>

namespace po = boost::program_options;
using namespace std;

int main(int argc, char **argv){
	int vpr = 10;
	string input_path;
	double x = 0;
	double y = 0;
	po::options_description desc("query usage");
	desc.add_options()
		("help,h", "produce help message")
		("partition,p", "use partition for distance calculating")
		("vertices_per_raster,n", po::value<int>(&vpr), "number of vertices per raster")
		("source_path,s", po::value<string>(&input_path)->required(), "path to the source polygons")
		("longitude,x", po::value<double>(&x), "longitude of query point")
		("latitude,y", po::value<double>(&y), "latitude of query point")

		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
		cout << desc << "\n";
		return 0;
	}
	po::notify(vm);
	timeval start = get_cur_time();
	query_context ctx;
	vector<MyPolygon *> polys = MyPolygon::load_binary_file(input_path.c_str(), ctx);
	logt("read polygons", start);
	for(MyPolygon *p:polys){
		p->partition(vpr);
	}
	logt("partitioning polygons", start);
	Point target(x,y);
	query_context ctx;
	ctx.use_grid = false;
	for(MyPolygon *p:polys){
		double dist = p->distance(target,&ctx);
		//printf("%f\n",dist);
	}
	logt("querying without rasterization", start);

	ctx.use_grid = true;
	for(MyPolygon *p:polys){
		p->distance(target,&ctx);
	}
	logt("querying with rasterization", start);

	for(MyPolygon *p:polys){
		delete p;
	}
	return 0;
}


