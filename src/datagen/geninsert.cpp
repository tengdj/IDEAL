#include "../geometry/MyPolygon.h"
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <vector>
#include <string.h>
#include <fstream>
#include <iostream>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

using namespace std;


int main(int argc, char** argv) {
	string inpath;
	int small_threshold = 500;
	int big_threshold = 1000000;
	po::options_description desc("load usage");
	string del("|");
	desc.add_options()
		("help,h", "produce help message")
		("path,p", po::value<string>(&inpath)->required(), "source path")
		("delimiter,d", po::value<string>(&del), "delimiter")
		("small_threshold,s", po::value<int>(&small_threshold), "minimum number of vertices for big polygons")
		("big_threshold,b", po::value<int>(&big_threshold), "maximum number of vertices for big polygons")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
		cout << desc << "\n";
		return 0;
	}
	po::notify(vm);

	ifstream infile;
	infile.open(inpath.c_str(), ios::in | ios::binary);
	int index = 0;
	while(!infile.eof()){
		MyPolygon * poly = MyPolygon::read_polygon_binary_file(infile);
		if(!poly||poly->get_num_vertices()<small_threshold||poly->get_num_vertices()>big_threshold){
			continue;
		}
		poly->setid(++index);
		printf("%d%s\"",poly->getid(),del.c_str());
		poly->print_without_return(false);
		printf("\"\n");
		delete poly;
	}
	infile.close();
}
