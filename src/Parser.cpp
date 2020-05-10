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
	//polygons[0]->print();
	mpoly->print();
	delete mpoly;
	return 0;
}


