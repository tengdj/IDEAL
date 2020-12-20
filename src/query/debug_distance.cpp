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
	Point p(-117.156449,33.141541);
	query_context ctx;
	ctx.use_grid = true;
	ctx.distance_buffer_size = 10;
	ctx.query_type = QueryType::within;
	polygon->distance(p, &ctx);
	delete polygon;
	return 0;
}




