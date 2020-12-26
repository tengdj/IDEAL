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

	query_context global_ctx;
	global_ctx = get_parameters(argc, argv);

	global_ctx.source_polygons = MyPolygon::load_binary_file(global_ctx.source_path.c_str(),global_ctx);
	if(global_ctx.use_grid||global_ctx.use_qtree){
		process_partition(&global_ctx);
	}
	if(global_ctx.use_convex_hull){
		process_convex_hull(&global_ctx);
	}

	if(global_ctx.use_mer){
		process_mer(&global_ctx);
	}


	size_t data_size = 0;
	size_t partition_size = 0;
	for(MyPolygon *poly:global_ctx.source_polygons){
		data_size += poly->get_data_size();
		if(global_ctx.use_grid){
			partition_size += poly->partition_size();
		}
		if(global_ctx.use_qtree){
			const int lfc = poly->get_qtree()->leaf_count();
			int tlfc = lfc;
			int bits = 1;
			while(tlfc>0){
				bits++;
				tlfc /= 2;
			}
			int bits_per_child = (bits+7)/8*8;
			partition_size += lfc*(bits_per_child*4+2)/8;
		}
		if(global_ctx.use_convex_hull){
			partition_size += poly->convex_hull->get_data_size();
		}
		if(global_ctx.use_mer){
			partition_size += 4*8;
		}
	}
	printf("%f\n",partition_size*100.0/data_size);

	return 0;
}


