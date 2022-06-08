#include "../include/MyPolygon.h"
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
	po::options_description desc("load usage");
	string del("|");
	string table_name;
	query_context ctx;
	desc.add_options()
		("help,h", "produce help message")
		("path,p", po::value<string>(&ctx.source_path)->required(), "source path")
		("delimiter,d", po::value<string>(&del), "delimiter")
		("table_name,t", po::value<string>(&table_name), "table name")
		("sql","generate insert sql script")
		("sample_rate,r", po::value<float>(&ctx.sample_rate), "max number of points inserted")
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

	int index = 0;
	if(!vm.count("point")){
		vector<MyPolygon *> polygons = MyPolygon::load_binary_file(ctx.source_path.c_str(), ctx, false);
		for(MyPolygon *poly:polygons){
			if(!tryluck(ctx.sample_rate)){
				if(poly){
					delete poly;
				}
				continue;
			}
			poly->setid(++index);

			if(vm.count("sql")){
				printf("INSERT INTO %s(id, geog) ",table_name.c_str());
				printf("VALUES (%ld, geography::STGeomFromText('%s', 4326).MakeValid());\n",poly->getid(), poly->to_string(true, true).c_str());
				if(index%100==0){
					printf("go\n");
				}
			}else{
				printf("%ld'%s",poly->getid(),del.c_str());
				poly->print_without_return(true, true);
				printf("'\n");
			}
			delete poly;
		}
	}else{
		ctx.target_path = ctx.source_path;
		ctx.load_points();
		for(int i=0;i<ctx.target_num;i++){
			if(!tryluck(ctx.sample_rate)){
				continue;
			}
			Point p(ctx.points[2*i],ctx.points[2*i+1]);
			index++;
			if(vm.count("sql")){
				printf("INSERT INTO %s(id, geog) ",table_name.c_str());
				printf("VALUES (%d, geography::STGeomFromText('%s', 4326));\n",index, p.to_string().c_str());
			}else{
				printf("%d'%s",index,del.c_str());
				p.print_without_return();
				printf("'\n");
			}
		}
	}

	log("generated %d",index);
}
