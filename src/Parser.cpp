/*
 * Parser.cpp
 *
 *  Created on: May 9, 2020
 *      Author: teng
 */

#include "MyPolygon.h"




int main(int argc, char **argv){
	int dim = 20;
	if(argc>1){
		dim = atoi(argv[1]);
	}

	string input_line;
	getline(cin, input_line);
	vector<string> result;
	tokenize(input_line, result, "\t");
	timeval start = get_cur_time();
	MyMultiPolygon *mpoly = new MyMultiPolygon(result[0].c_str());
	logt("parse polygon", start);
	MyPolygon *poly = mpoly->get_polygon(0);

	vector<vector<Pixel>> partitions = poly->partition(dim, dim);
	logt("partitioning polygon", start);

	MyMultiPolygon *inpolys = new MyMultiPolygon();
	MyMultiPolygon *borderpolys = new MyMultiPolygon();
	MyMultiPolygon *outpolys = new MyMultiPolygon();
	char *data = MyPolygon::encode(partitions);
	logt("encoding", start);
	vector<vector<Pixel>> newpartitions = MyPolygon::decode(data);
	logt("decoding", start);

	for(int i=0;i<dim;i++){
		for(int j=0;j<dim;j++){
			MyPolygon *m = newpartitions[i][j].to_polygon();
			if(newpartitions[i][j].status==BORDER){
				borderpolys->insert_polygon(m);
			}else if(newpartitions[i][j].status==IN){
				inpolys->insert_polygon(m);
			}else if(newpartitions[i][j].status==OUT){
				outpolys->insert_polygon(m);
			}
		}
	}
	logt("allocating partitions", start);
	Point p(-120.7985,38.6137);
	for(int i=0;i<10000;i++){
		poly->contain(p,false);
	}
	logt("querying", start);

//	poly->print();
//	cout<<"border"<<endl;
//	borderpolys->print();
//	cout<<"in"<<endl;
//	inpolys->print();
//	cout<<"out"<<endl;
//	outpolys->print();
	delete mpoly;
	delete borderpolys;
	delete inpolys;
	delete outpolys;
	return 0;
}


