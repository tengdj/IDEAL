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
	int max_dimx = 40;
	int max_dimy = 40;
	po::options_description desc("query usage");
	desc.add_options()
		("help,h", "produce help message")
		("max_dimx,x", po::value<int>(&max_dimx), "max dimension on horizontal")
		("max_dimy,y", po::value<int>(&max_dimy), "max dimension on vertical")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
		cout << desc << "\n";
		return 0;
	}
	po::notify(vm);
	timeval start = get_cur_time();
	MyPolygon *poly = MyMultiPolygon::read_one_polygon();
	MyPolygon *poly2 = poly->clone();
	logt("read polygon", start);
	cout<<poly->get_num_vertices()<<endl;
	vector<vector<Pixel>> partitions = poly->partition(max_dimx, max_dimy);
	logt("partitioning polygon", start);
	vector<vector<Pixel>> partitions2 = poly2->partition_with_query(max_dimx, max_dimy);
	logt("partitioning polygon with query", start);
//	assert(partitions.size()==partitions2.size());
//	for(int i=0;i<partitions.size();i++){
//		assert(partitions[i].size()==partitions2[i].size());
//		for(int j=0;j<partitions2.size();j++){
//			assert(partitions[i][j].status==partitions2[i][j].status);
//		}
//	}


	delete poly;
	return 0;
}


