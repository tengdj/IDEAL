/*
 * debug_distance.cpp
 *
 *  Created on: Dec 19, 2020
 *      Author: teng
 */

#include "../geometry/MyPolygon.h"
#include <fstream>
#include <boost/program_options.hpp>
#include <queue>

namespace po = boost::program_options;
using namespace std;

int main(int argc, char **argv){


	query_context ctx;
	Point p(-117.156449,33.141541);

//	MyPolygon *polygon = MyPolygon::load_binary_file_single("/data/gisdata/raster/dat/all_source_100.dat", ctx, 198009);
	MyPolygon *polygon=MyMultiPolygon::read_one_polygon();
	polygon->triangulate();

	polygon->print_triangles();

//	polygon->print();
//	return 0;
//	polygon->get_convex_hull();
//	polygon->partition(10);
//	ctx.use_triangulate = true;
//	polygon->valid_for_triangulate = true;
//	polygon->triangulate();
//	polygon->build_rtree();

//	p.print();
//	polygon->print();
//	polygon->print_triangles();

//	polygon->build_rtree();
//	Pixel *root = polygon->get_rtree_pixel();
//	queue<Pixel *> pixq;
//	pixq.push(root);
//	int level = 0;
//	while(pixq.size()>0){
//		int s = pixq.size();
//		printf("level %d %d\n",level++,s);
//		printf("MULTIPOLYGON(");
//		for(int i=0;i<s;i++){
//			Pixel *cur = pixq.front();
//			pixq.pop();
//			if(i>0){
//				printf(",");
//			}
//			printf("((");
//			cur->print_vertices();
//			printf("))");
//			for(Pixel *p:cur->children){
//				pixq.push(p);
//			}
//		}
//		printf(")\n");
//	}

	double d = polygon->distance(p, &ctx);
	printf("%f\n",d);

//	ctx.use_grid = true;
//	ctx.distance_buffer_size = 10;
//	ctx.query_type = QueryType::within;
//	ctx.mer_sample_round = 20;
//	polygon->distance(p, &ctx);
//
//	VertexSequence *ch = polygon->get_convex_hull();
//
//
//
//	polygon->print();
//	query_context qt;
//	qt.use_grid = true;
//	polygon->print_partition(qt);
//
//	assert(ch);
//	cout<<"Polygon (";
//	ch->print();
//	cout<<")"<<endl;
//
//	Pixel *mbr = polygon->getMBB();
//
//	Pixel *pix = polygon->getMER(&ctx);
//	MyPolygon::gen_box(*pix)->print();
//
	delete polygon;
	return 0;
}




