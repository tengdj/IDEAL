/*
 * MyPolygon.cpp
 *
 *  Created on: May 9, 2020
 *      Author: teng
 */

#include "MyPolygon.h"
#include <string.h>
#include <assert.h>
#include <iostream>


inline bool is_number(char ch){
	return ch=='-'||ch=='.'||(ch<='9'&&ch>='0');
}

inline double read_double(const char *input, size_t &offset){
	char tmp[100];
	while(!is_number(input[offset])){
		offset++;
	}
	int index = 0;
	while(is_number(input[offset])){
		tmp[index++] = input[offset++];
	}
	tmp[index] = '\0';
	return atof(tmp);
}
inline void skip_space(const char *input, size_t &offset){
	while(input[offset]==' '||input[offset]=='\t'||input[offset]=='\n'){
		offset++;
	}
}


double *MyPolygon::read_vertices(const char *wkt, size_t &offset){
	// read until the left parenthesis
	skip_space(wkt, offset);
	assert(wkt[offset++]=='(');

	// count the number of vertices
	int cur_offset = offset;
	int num_vertices = 0;
	while(wkt[cur_offset++]!=')'){
		if(wkt[cur_offset]==','){
			num_vertices++;
		}
	}
	num_vertices++;
	double *vertices = new double[2*num_vertices+1];
	vertices[0] = num_vertices;
	//
	int index = 1;
	for(int i=0;i<num_vertices*2;i++){
		vertices[index++] = read_double(wkt, offset);
	}
	// read until the right parenthesis
	skip_space(wkt, offset);
	assert(wkt[offset++]==')');
	return vertices;
}

inline double *copy_vertices(double *from){
	assert(from);
	double num_vertices = *from;
	double *boundary = new double[(int)num_vertices*2+1];
	memcpy((void *)boundary,(void *)from,sizeof(double)*((int)num_vertices*2+1));
	return boundary;
}

MyPolygon *MyPolygon::clone(){
	MyPolygon *polygon = new MyPolygon();
	polygon->boundary = copy_vertices(boundary);
	for(double *t:internal_polygons){
		polygon->internal_polygons.push_back(copy_vertices(t));
	}
	return polygon;
}


MyPolygon *MyPolygon::read_polygon(const char *wkt, size_t &offset){
	MyPolygon *polygon = new MyPolygon();
	skip_space(wkt, offset);
	// left parentheses for the entire polygon
	assert(wkt[offset++]=='(');

	// read the vertices of the boundary polygon
	polygon->boundary = read_vertices(wkt, offset);

	skip_space(wkt, offset);
	//polygons as the holes of the boundary polygon
	while(wkt[offset]==','){
		offset++;
		polygon->internal_polygons.push_back(read_vertices(wkt, offset));
		skip_space(wkt, offset);
	}
	assert(wkt[offset++]==')');
	return polygon;
}

MyPolygon::~MyPolygon(){
	if(this->boundary){
		delete []boundary;
	}
	for(double *p:internal_polygons){
		if(p){
			delete []p;
		}
	}
	internal_polygons.clear();
}

void MyPolygon::print_without_head(){
	cout<<"((";
	double *cur_b = boundary;
	cur_b++;
	for(int i=0;i<boundary[0];i++){
		if(i!=0){
			cout<<",";
		}
		printf("%f ",*cur_b);
		cur_b++;
		printf("%f",*cur_b);
		cur_b++;
	}
	cout<<")";
	for(double *iboundary:this->internal_polygons){
		cout<<", (";
		cur_b = iboundary;
		cur_b++;
		for(int i=0;i<iboundary[0];i++){
			if(i!=0){
				cout<<",";
			}
			printf("%f ",*cur_b);
			cur_b++;
			printf("%f",*cur_b);
			cur_b++;
		}
		cout<<")";
	}
	cout<<")";
}

void MyPolygon::print(){
	cout<<"POLYGON";
	print_without_head();
	cout<<endl;
}

MyMultiPolygon::MyMultiPolygon(const char *wkt){
	size_t offset = 0;
	// read the symbol MULTIPOLYGON
	while(wkt[offset]!='M'){
		offset++;
	}
	for(int i=0;i<strlen(multipolygon_char);i++){
		assert(wkt[offset++]==multipolygon_char[i]);
	}
	skip_space(wkt,offset);
	// read the left parenthesis
	assert(wkt[offset++]=='(');

	// read the first polygon
	polygons.push_back(MyPolygon::read_polygon(wkt,offset));
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



