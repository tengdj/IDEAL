/*
 * Parser.cpp
 *
 *  Created on: May 9, 2020
 *      Author: teng
 */

#include "../include/MyPolygon.h"
#include <fstream>
#include "../index/RTree.h"
#include <queue>
#include <boost/program_options.hpp>

namespace po = boost::program_options;
using namespace std;

RTree<MyPolygon *, double, 2, double> tree;

bool MySearchCallback(MyPolygon *poly, void* arg){
	query_context *ctx = (query_context *)arg;
	Point p = *(Point *)ctx->target;
	// MBB check
	if(!poly->getMBB()->contain(p)){
		return true;
	}
	if(poly->contain(p, ctx)){
		poly->print();
	}

	// keep going until all hit objects are found
	return true;
}


int main(int argc, char** argv) {
	string source_path;
	double point[2];
	po::options_description desc("query usage");
	desc.add_options()
		("help,h", "produce help message")
		("source,s", po::value<string>(&source_path)->required(), "path to the source")
		("x_value,x", po::value<double>(&point[0])->required(), "x value for target")
		("y_value,y", po::value<double>(&point[1])->required(), "y value for target")
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
	vector<MyPolygon *> source = load_binary_file(source_path.c_str(),ctx);
	logt("loaded %ld polygons", start, source.size());
	for(MyPolygon *p:source){
		p->rasterization(10);
	}
	logt("partition polygons", start);

	for(MyPolygon *p:source){
		tree.Insert(p->getMBB()->low, p->getMBB()->high, p);
	}
	logt("building R-Tree", start);
	ctx.use_grid = true;
	Point p(point[0],point[1]);
	ctx.target = (void *)&p;
	tree.Search(point, point, MySearchCallback, (void *)&ctx);

	for(MyPolygon *p:source){
		delete p;
	}

	return 0;
}



