/*
 * Parser.cpp
 *
 *  Created on: May 9, 2020
 *      Author: teng
 */

#include "../geometry/MyPolygon.h"
#include <fstream>
#include <boost/program_options.hpp>
#include <stack>

namespace po = boost::program_options;
using namespace std;

int main(int argc, char **argv){

	po::options_description desc("query usage");
	query_context ctx;

	desc.add_options()
		("help,h", "produce help message")
		("print,p", "print the pixels")
		("quad_tree,q", "partition with Q-tree")
		("vertices_per_raster,v", po::value<int>(&ctx.vpr), "number of vertices per raster")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
		cout << desc << "\n";
		return 0;
	}
	po::notify(vm);
	ctx.use_qtree = vm.count("quad_tree");
	ctx.use_grid = !ctx.use_qtree;

	MyPolygon *poly = MyMultiPolygon::read_one_polygon();
	if(ctx.use_grid){
		poly->partition(ctx.vpr);
	}
	if(ctx.use_qtree){
		poly->partition_qtree(ctx.vpr);
	}
	if(vm.count("print")){
		poly->print_partition(ctx);
		ctx.use_grid = true;
		ctx.use_qtree = false;
		//poly->print_partition(ctx);
	}


	delete poly;
	return 0;
}


