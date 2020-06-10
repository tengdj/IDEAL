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
	string path;

	po::options_description desc("load usage");
	desc.add_options()
		("help,h", "produce help message")
		("output,o", po::value<string>(&path)->required(), "path to the output file")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
		cout << desc << "\n";
		return 0;
	}
	po::notify(vm);
	ofstream os;
	os.open(path.c_str(), ios::out | ios::binary |ios::trunc);

	struct timeval start_time = get_cur_time();
	std::string input_line;
	int num_objects = 0;
	const int buffer_size = 10000;
	double *buffer = new double[buffer_size*2*sizeof(double)];
	int index = 0;

	while(getline(cin, input_line)) {
		Point *p = Point::read_one_point(input_line);
		if(!p){
			continue;
		}
		buffer[index++] = p->x;
		buffer[index++] = p->y;
		if(index==buffer_size*2){
			os.write((char *)buffer, index*sizeof(double));
			index = 0;
		}
		num_objects++;
	}
	if(index>0){
		os.write((char *)buffer, index*sizeof(double));
	}
	logt("processed %d objects", start_time, num_objects);
	os.close();
	delete []buffer;
}
