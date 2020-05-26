/*
 * MyMultipolygon.cpp
 *
 *  Created on: May 26, 2020
 *      Author: teng
 */

#include "MyPolygon.h"

MyMultiPolygon::MyMultiPolygon(const char *wkt){
	size_t offset = 0;
	// read the symbol MULTIPOLYGON
	while(wkt[offset]!='M'&&wkt[offset]!='P'){
		offset++;
	}
	bool is_multiple = (wkt[offset]=='M');
	if(is_multiple){
		for(int i=0;i<strlen(multipolygon_char);i++){
			assert(wkt[offset++]==multipolygon_char[i]);
		}
		skip_space(wkt,offset);
		// read the left parenthesis
		assert(wkt[offset++]=='(');
	}else {
		for(int i=0;i<strlen(polygon_char);i++){
			assert(wkt[offset++]==polygon_char[i]);
		}
		skip_space(wkt,offset);
	}
	// read the first polygon
	polygons.push_back(MyPolygon::read_polygon(wkt,offset));
	if(is_multiple){
		skip_space(wkt, offset);
		// read the rest polygons
		while(wkt[offset]==','){
			offset++;
			MyPolygon *in = MyPolygon::read_polygon(wkt,offset);
			assert(in);
			polygons.push_back(in);
			skip_space(wkt,offset);
		}
		// read the right parenthesis
		assert(wkt[offset++]==')');
	}
}

MyMultiPolygon::~MyMultiPolygon(){
	for(MyPolygon *p:polygons){
		delete p;
	}
}

void MyMultiPolygon::print(){
	cout<<"MULTIPOLYGON (";
	for(int i=0;i<polygons.size();i++){
		if(i>0){
			cout<<",";
		}
		polygons[i]->print_without_head();
	}
	cout<<")"<<endl;
}

void MyMultiPolygon::print_separate(){
	for(MyPolygon *p:polygons){
		p->print();
	}
}



