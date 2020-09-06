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
	int small_threshold = 500;
	int big_threshold = 1000000;
	po::options_description desc("load usage");
	string del("|");
	string table_name;
	query_context qc;
	int max_count = -1;
	desc.add_options()
		("help,h", "produce help message")
		("path,p", po::value<string>(&qc.source_path)->required(), "source path")
		("delimiter,d", po::value<string>(&del), "delimiter")
		("table_name,t", po::value<string>(&table_name), "table name")
		("sql","generate insert sql script")
		("small_threshold,s", po::value<int>(&small_threshold), "minimum number of vertices for big polygons")
		("big_threshold,b", po::value<int>(&big_threshold), "maximum number of vertices for big polygons")
		("max_number,m", po::value<int>(&max_count), "max number of points inserted")
		("point","the source is point")

		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
		cout << desc << "\n";
		return 0;
	}
	if(vm.count("sql")&&!vm.count("table_name")){
		cout<<"table_name should be specified for sql generation"<<endl;
		cout << desc << "\n";
		return 0;
	}
	po::notify(vm);

	if(vm.count("sql")){
		printf("SET QUOTED_IDENTIFIER ON;\n");
	}

	if(!vm.count("point")){
		ifstream infile;
		infile.open(qc.source_path.c_str(), ios::in | ios::binary);
		int index = 0;
		while(!infile.eof()){
			MyPolygon * poly = MyPolygon::read_polygon_binary_file(infile);
			if(!poly||poly->get_num_vertices()<small_threshold||poly->get_num_vertices()>big_threshold){
				continue;
			}
			poly->setid(++index);

			if(vm.count("sql")){
				printf("INSERT INTO %s(id, geog) ",table_name.c_str());
				printf("VALUES (%d, geography::STGeomFromText('%s', 4326).MakeValid());\n",poly->getid(), poly->to_string().c_str());
				if((index)%1000==0){
					printf("go;\n");
				}
			}else{
				printf("%d%s'",poly->getid(),del.c_str());
				poly->print_without_return(false);
				printf("'\n");
			}
			delete poly;
			if(max_count>0&&index==max_count){
				break;
			}
		}
		infile.close();
	}else{
		qc.target_path = qc.source_path;
		qc.load_points();
		if(max_count==-1||max_count>qc.target_num){
			max_count = qc.target_num;
		}
		for(int i=0;i<max_count;i++){
			Point p(qc.points[2*i],qc.points[2*i+1]);
			if(vm.count("sql")){
				printf("INSERT INTO %s(id, geog) ",table_name.c_str());
				printf("VALUES (%d, geography::STGeomFromText('%s', 4326).MakeValid());\n",i+1, p.to_string().c_str());
				if((i+1)%1000==0){
					printf("go;\n");
				}
			}else{
				printf("%d%s'",i+1,del.c_str());
				p.print_without_return();
				printf("'\n");
			}

		}
	}
}
