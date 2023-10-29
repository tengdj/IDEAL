/*
 * Parser.cpp
 *
 *  Created on: May 9, 2020
 *      Author: teng
 */



#include "../index/RTree.h"
#include <queue>
#include <fstream>
#include "../include/MyPolygon.h"
#include <omp.h>
// some shared parameters

int main(int argc, char** argv) {

	query_context cctx;
	cctx.query_type = QueryType::distance;
	cctx.geography = false;
	vector<MyPolygon *> polys = load_binary_file(argv[1], cctx);
	int vpr = atoi(argv[2]);

	vector<vector<double>> values;
	size_t total = 0;
	for(int i=0;i<20;i++){
		vector<double> v;
		values.push_back(v);
	}
#pragma omp parallel for
	for(int i=0;i<polys.size();i++){
		int num = polys[i]->get_num_vertices();
		if(num<1000 || num>=20000){
			continue;
		}
		polys[i]->rasterization(vpr);
		double it = polys[i]->get_rastor()->get_num_intersection();

//#pragma omp critical
//		if(it<1){
//			polys[i]->print();
//			polys[i]->get_rastor()->print();
//			cout<<polys[i]->get_rastor()->get_dimx()<<endl;
//			cout<<it<<endl;
//
//			exit(0);
//		}

		if(num>=20000){
			continue;
		}
		num/=1000;
#pragma omp critical
		{
			total++;
			values[num].push_back(it);
		}
	}
	for(int i=0;i<20;i++){
		if(values[i].size()==0){
			continue;
		}
		sort(values[i].begin(), values[i].end());
		double avg = 0;
		int lst = 0;
		int fst = values[i].size()*0.01;
		size_t total = 0;
		for(int j=fst;j<values[i].size()*0.99;j++){
			avg += values[i][j];
			total++;
			lst = j;
		}
		printf("%d\t%ld\t%f\t%f\t%f\t%f\n",i*1000,total,values[i][0],values[i][lst/2],avg/total,values[i][lst]);
	}


	return 0;
}



