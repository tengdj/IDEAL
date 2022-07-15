/*
 * Parser.cpp
 *
 *  Created on: May 9, 2020
 *      Author: teng
 */

#include <fstream>
#include <boost/program_options.hpp>

#include "MyPolygon.h"
#include "mygpu.h"

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
		("gpu,g", "run with gpu")
		("vertices_per_raster,v", po::value<int>(&vpr), "number of vertices per raster")
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
	ctx.gpu = vm.count("gpu");
	vector<MyPolygon *> polys = load_binary_file(input_path.c_str(), ctx);
	logt("read polygons", start);

	Point target(x,y);

	init_gpus();

	double tdist = 0;
	start = get_cur_time();
	ctx.gpu = false;
	for(MyPolygon *p:polys){
		double dist = p->distance(target,&ctx);
		tdist += dist;
		//logt("%f\t%d",start, dist, p->get_num_vertices());
		//break;
	}
	logt("querying without gpu %f", start,tdist);
	ctx.gpu = true;
	for(MyPolygon *p:polys){
		double dist = p->distance(target,&ctx);
		tdist += dist;
		//logt("%f\t%d",start, dist, p->get_num_vertices());
		//break;
	}
	logt("querying with gpu %f", start, tdist);

	if(vm.count("vertices_per_raster")){
		for(MyPolygon *p:polys){
			p->rasterization(vpr);
		}
		logt("rasterizing polygons", start);

		ctx.use_grid = true;
		for(MyPolygon *p:polys){
			p->distance(target,&ctx);
		}
		logt("querying with rasterization", start);
	}

	for(MyPolygon *p:polys){
		delete p;
	}
	release_gpus();
	return 0;
}


