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
	int vpr = 10;
	po::options_description desc("query usage");
	desc.add_options()
		("help,h", "produce help message")
		("print,p", "print the pixels")
		("quad_tree,q", "partition with Q-tree")
		("scanline,s", "partition with scanline")
		("vertices_per_raster,v", po::value<int>(&vpr), "number of vertices per raster")
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
	vector<vector<Pixel>> partitions;
	if(vm.count("scanline")){
		partitions = poly->partition_scanline(vpr);
		logt("partitioning polygon with scanline", start);
	}else if(vm.count("quad_tree")){
		poly->partition_qtree(vpr);
		logt("partitioning polygon with Q-tree", start);
	}else{
		partitions = poly->partition(vpr);
		logt("partitioning polygon %d %d", start,partitions.size(),partitions[0].size());
	}

	MyMultiPolygon *inpolys = new MyMultiPolygon();
	MyMultiPolygon *borderpolys = new MyMultiPolygon();
	MyMultiPolygon *outpolys = new MyMultiPolygon();

	if(vm.count("quad_tree")){
		QTNode *qtree = poly->get_qtree();
		log("leaf count %d",qtree->leaf_count());
		log("size in bytes %d",qtree->size());
		std::stack<QTNode *> ws;
		ws.push(qtree);
		while(!ws.empty()){
			QTNode *cur = ws.top();
			ws.pop();
			if(cur->isleaf){
				if(cur->interior){
					MyPolygon *m = cur->mbb.to_polygon();
					inpolys->insert_polygon(m);
				}else if(cur->exterior){
					MyPolygon *m = cur->mbb.to_polygon();
					outpolys->insert_polygon(m);
				}else{
					MyPolygon *m = cur->mbb.to_polygon();
					borderpolys->insert_polygon(m);
				}
			}else{
				cur->push(ws);
			}
		}
	}else{
		for(int i=0;i<partitions.size();i++){
			for(int j=0;j<partitions[0].size();j++){
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


