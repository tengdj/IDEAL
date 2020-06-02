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
		("print,p", "print the pixels")
		("query,q", "partition with query")
		("scanline,s", "partition with scanline")
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
	logt("read one polygon", start);
	vector<vector<Pixel>> partitions = poly->partition(max_dimx, max_dimy);
	logt("partitioning polygon", start);
	if(vm.count("scanline")){
		poly->reset_partition();
		partitions = poly->partition_scanline(max_dimx, max_dimy);
		logt("partitioning polygon with scanline", start);
	}
	if(vm.count("query")){
		poly->reset_partition();
		partitions = poly->partition_with_query(max_dimx, max_dimy);
		logt("partitioning polygon with query", start);
	}

	MyMultiPolygon *inpolys = new MyMultiPolygon();
	MyMultiPolygon *borderpolys = new MyMultiPolygon();
	MyMultiPolygon *outpolys = new MyMultiPolygon();

	for(int i=0;i<max_dimx;i++){
		for(int j=0;j<max_dimy;j++){
			MyPolygon *m = partitions[i][j].to_polygon();
			if(partitions[i][j].status==BORDER){
				borderpolys->insert_polygon(m);
			}else if(partitions[i][j].status==IN){
				inpolys->insert_polygon(m);
			}else if(partitions[i][j].status==OUT){
				outpolys->insert_polygon(m);
			}
		}
	}
	logt("allocating partitions", start);

	if(vm.count("print")){
		cout<<"polygon:"<<endl;
		poly->print();
		cout<<"border:"<<endl;
		borderpolys->print();
		cout<<"in:"<<endl;
		inpolys->print();
		cout<<"out:"<<endl;
		outpolys->print();
	}

	delete borderpolys;
	delete inpolys;
	delete outpolys;
	delete poly;
	return 0;
}


