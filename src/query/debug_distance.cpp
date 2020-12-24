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
	ctx.mer_sample_round = 20;
	polygon->distance(p, &ctx);

	VertexSequence *ch = polygon->get_convex_hull();



	polygon->print();
	query_context qt;
	qt.use_grid = true;
	polygon->print_partition(qt);

	assert(ch);
	cout<<"Polygon (";
	ch->print();
	cout<<")"<<endl;

	Pixel *mbr = polygon->getMBB();

	Pixel *pix = polygon->getMER(&ctx);
	pix->to_polygon()->print();

	delete polygon;
	return 0;
}




