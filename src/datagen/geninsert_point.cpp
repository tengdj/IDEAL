#include "../geometry/MyPolygon.h"
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <vector>
#include <string.h>
#include <fstream>
#include <iostream>
using namespace std;


int main(int argc, char** argv) {

	int max_count = -1;
	po::options_description desc("load usage");
	string del("|");
	query_context qc;
	desc.add_options()
		("help,h", "produce help message")
		("path,p", po::value<string>(&qc.target_path)->required(), "source path")
		("delimiter,d", po::value<string>(&del), "delimiter")
		("max_number,m", po::value<int>(&max_count), "max number of points inserted")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
		cout << desc << "\n";
		return 0;
	}
	po::notify(vm);

	qc.load_points();
	if(max_count==-1||max_count>qc.target_num){
		max_count = qc.target_num;
	}
	for(int i=0;i<max_count;i++){
		Point p(qc.points[2*i],qc.points[2*i+1]) ;
		printf("%d%s\"",i+1,del.c_str());
		p.print_without_return();
		printf("\"\n");
	}
}
