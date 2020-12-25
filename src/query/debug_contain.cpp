/*
 * debug_distance.cpp
 *
 *  Created on: Dec 19, 2020
 *      Author: teng
 */

#include "../geometry/MyPolygon.h"
#include <fstream>
#include <boost/program_options.hpp>

namespace po = boost::program_options;
using namespace std;

int main(int argc, char **argv){
	MyPolygon *polygon=MyMultiPolygon::read_one_polygon();
	polygon->partition(10);
	Point p(-117.915111,33.799661);
	query_context ctx;
	ctx.use_grid = true;
	ctx.query_type = QueryType::contain;


	polygon->print();
	query_context qt;
	qt.use_grid = true;
	polygon->print_partition(qt);

	cout<<polygon->contain(p, &ctx)<<endl;
	delete polygon;
	return 0;
}




