/*
 * Parser.cpp
 *
 *  Created on: May 9, 2020
 *      Author: teng
 */

#include "MyPolygon.h"




int main(){
	string input_line;
	getline(cin, input_line);
	vector<string> result;
	tokenize(input_line, result, "\t");
	MyMultiPolygon *mpoly = new MyMultiPolygon(result[0].c_str());
	MyPolygon *poly = mpoly->get_polygon(0);
	MyPolygon *mbb = poly->getMBB();
	vector<MyPolygon *> parts = poly->partition(10, 10);
	MyMultiPolygon *mypolys = new MyMultiPolygon();
	mypolys->insert_polygon(poly->clone());
	mypolys->insert_polygon(parts);
	mypolys->print();
	delete mpoly;
	delete mypolys;
	return 0;
}


