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
#include <float.h>


inline bool is_number(char ch){
	return ch=='-'||ch=='.'||(ch<='9'&&ch>='0')||ch=='e';
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
	if(mbb){
		delete mbb;
	}
}


void Pixel::process_enter_leave(Direction enter_d, double enter_val, Direction leave_d, double leave_val){
	// totally 16 possible combinations
	if(enter_d==LEFT&&leave_d==LEFT){
		left = true;
		if(leave_val>enter_val){
			right = true;
			top = true;
			bottom = true;
		}
	}else if(enter_d==LEFT&&leave_d==TOP){
		left = true;
		right = true;
		top = true;
		bottom = true;
	}else if(enter_d==LEFT&&leave_d==RIGHT){
		left = true;
		right = true;
		bottom = true;
	}else if(enter_d==LEFT&&leave_d==BOTTOM){
		left = true;
		bottom = true;
	}else if(enter_d==TOP&&leave_d==LEFT){
		left = true;
		top = true;
	}else if(enter_d==TOP&&leave_d==TOP){
		top = true;
		if(leave_val>enter_val){
			left = true;
			right = true;
			bottom = true;
		}
	}else if(enter_d==TOP&&leave_d==RIGHT){
		left = true;
		right = true;
		top = true;
		bottom = true;
	}else if(enter_d==TOP&&leave_d==BOTTOM){
		left = true;
		top = true;
		bottom = true;
	}else if(enter_d==RIGHT&&leave_d==LEFT){
		left = true;
		right = true;
		top = true;
	}else if(enter_d==RIGHT&&leave_d==TOP){
		right = true;
		top = true;
	}else if(enter_d==RIGHT&&leave_d==RIGHT){
		right = true;
		if(leave_val<enter_val){
			left = true;
			top = true;
			bottom = true;
		}
	}else if(enter_d==RIGHT&&leave_d==BOTTOM){
		left = true;
		right = true;
		top = true;
		bottom = true;
	}else if(enter_d==BOTTOM&&leave_d==LEFT){
		left = true;
		right = true;
		top = true;
		bottom = true;
	}else if(enter_d==BOTTOM&&leave_d==TOP){
		right = true;
		top = true;
		bottom = true;
	}else if(enter_d==BOTTOM&&leave_d==RIGHT){
		right = true;
		bottom = true;
	}else if(enter_d==BOTTOM&&leave_d==BOTTOM){
		bottom = true;
		if(leave_val<enter_val){
			left = true;
			right = true;
			top = true;
		}
	}else{
		assert(false);
	}
}

void Pixel::enter(double val, Direction d){
	if(true){
		cout<<direction_str[d];
		cout<<" enter "<<id[0]<<" "<<id[1]<<endl;
	}
	// one Pixel cannot be entered continuously twice
	assert(!entered);
	entered = true;
	if(leaved){
		process_enter_leave(d,val,in_direction, in_vertex);
	}else{
		in_vertex = val;
		in_direction = d;
	}
	this->status = BORDER;
}



void Pixel::leave(double val, Direction d){
	if(true){
		cout<<direction_str[d];
		cout<<" leave "<<id[0]<<" "<<id[1]<<endl;
	}

	//one Pixel cannot be left continuously twice
	assert(!leaved);
	leaved = true;
	if(entered){
		process_enter_leave(in_direction, in_vertex, d, val);
	}else{
		in_vertex = val;
		in_direction = d;
	}
	this->status = BORDER;
}

vector<MyPolygon *> MyPolygon::partition(const int dimx, const int dimy){

	assert(dimx>0);
	assert(dimy>0);
	vector<MyPolygon *> parts;
	getMBB();
	const double start_x = mbb->boundary[1];
	const double start_y = mbb->boundary[2];
	const double end_x = mbb->boundary[5];
	const double end_y = mbb->boundary[6];
	const double step_x = (end_x-start_x)/dimx;
	const double step_y = (end_y-start_y)/dimy;
	vector<vector<Pixel>> partitions;
	for(double i=0;i<=dimx;i++){
		vector<Pixel> v;
		for(double j=0;j<=dimy;j++){
			Pixel m;
			m.id[0] = i;
			m.id[1] = j;
			m.low[0] = i*step_x+start_x;
			m.high[0] = (i+1.0)*step_x+start_x;
			m.low[1] = j*step_y+start_y;
			m.high[1] = (j+1.0)*step_y+start_y;
			v.push_back(m);
		}
		partitions.push_back(v);
	}
	for(int i=0;i<num_boundary_vertices()-1;i++){

		double x1 = get_vertex_x(i);
		double y1 = get_vertex_y(i);
		double x2 = get_vertex_x(i+1);
		double y2 = get_vertex_y(i+1);

		int cur_startx = (x1-start_x)/step_x;
		int cur_endx = (x2-start_x)/step_x;
		int cur_starty = (y1-start_y)/step_y;
		int cur_endy = (y2-start_y)/step_y;

		if(cur_startx==dimx){
			cur_startx--;
		}
		if(cur_endx==dimx){
			cur_endx--;
		}

		int minx = min(cur_startx,cur_endx);
		int maxx = max(cur_startx,cur_endx);


		if(cur_starty==dimy){
			cur_starty--;
		}
		if(cur_endy==dimy){
			cur_endy--;
		}
		int miny = min(cur_starty,cur_endy);
		int maxy = max(cur_starty,cur_endy);

		partitions[cur_startx][cur_starty].status = BORDER;
		partitions[cur_endx][cur_endy].status = BORDER;

//		if(cur_endx==0&&cur_endy==2){
//			printf("%f %f %f %f\n",y1,cur_startx*step_x+start_x,x1,(x1-start_x)/step_x);
//			printf("%f %f %f %f\n",y2, cur_endx*step_x+start_x,x2,(x2-start_x)/step_x);
//		}


		//in the same pixel
		if(cur_startx==cur_endx&&cur_starty==cur_endy){
			continue;
		}

		if(y1==y2){
			//left to right
			if(cur_startx<cur_endx){
				for(int x=cur_startx;x<cur_endx;x++){
					partitions[x][cur_starty].leave(y1,RIGHT);
					partitions[x+1][cur_starty].enter(y1,LEFT);
				}
			}else { // right to left
				for(int x=cur_startx;x>cur_endx;x--){
					partitions[x][cur_starty].leave(y1, LEFT);
					partitions[x-1][cur_starty].enter(y1, RIGHT);
				}
			}
		}else if(x1==x2){
			//bottom up
			if(cur_starty<cur_endy){
				for(int y=cur_starty;y<cur_endy;y++){
					partitions[cur_startx][y].leave(x1, TOP);
					partitions[cur_startx][y+1].enter(x1, BOTTOM);
				}
			}else { //top down
				for(int y=cur_starty;y>cur_endy;y--){
					partitions[cur_startx][y].leave(x1, BOTTOM);
					partitions[cur_startx][y-1].enter(x1, TOP);
				}
			}
		}else{

			// todo: logic is wrong
			double a = (y1-y2)/(x1-x2);
			double b = (x1*y2-x2*y1)/(x1-x2);
			// cross with vertical lines
			if(cur_startx<cur_endx){//left to right
				for(int x=cur_startx;x<cur_endx;x++){
					double xval = ((double)x+1)*step_x+start_x;
					double yval = xval*a+b;
					int y = (yval-start_y)/step_y;
					partitions[x][y].leave(yval,RIGHT);
					partitions[x+1][y].enter(yval,LEFT);
				}
			}else{//right to left

				for(int x=cur_startx;x>cur_endx;x--){
					double xval = (double)x*step_x+start_x;
					double yval = xval*a+b;
					int y = (yval-start_y)/step_y;
					partitions[x][y].leave(yval,LEFT);
					partitions[x-1][y].enter(yval, RIGHT);
				}
			}

			// cross with horizontal lines
			if(cur_starty<cur_endy){//bottom up
				for(int y=cur_starty;y<cur_endy;y++){
					float yval = ((double)y+1)*step_y+start_y;
					float xval = (yval-b)/a;
					int x = (xval-start_x)/step_x;
					partitions[x][y].leave(xval, TOP);
					partitions[x][y+1].enter(xval, BOTTOM);
				}
			}else{//top down
				for(int y=cur_starty;y>cur_endy;y--){
					float yval = (double)y*step_y+start_y;
					float xval = (yval-b)/a;
					int x = (xval-start_x)/step_x;
					partitions[x][y].leave(xval, BOTTOM);
					partitions[x][y-1].enter(xval, TOP);
				}
			}
		}
	}

	//spread the in/out with border information
	for(int i=0;i<dimx;i++){
		for(int j=0;j<dimy;j++){
			// bottom up
			if(partitions[i][j].top){
				partitions[i][j+1].bottom = true;
				if(partitions[i][j+1].status==OUT){
					//partitions[i][j+1].setin();
				}
			}else if(partitions[i][j+1].bottom){
				partitions[i][j].top = true;
				if(partitions[i][j].status==OUT){
					//partitions[i][j].setin();
				}
			}
			// left to right
			if(partitions[i][j].right){
				partitions[i][j+1].left = true;
				if(partitions[i][j+1].status==OUT){
					//partitions[i][j+1].setin();
				}
			}else if(partitions[i][j+1].left){
				partitions[i][j].right = true;
				if(partitions[i][j].status==OUT){
					//partitions[i][j].setin();
				}
			}
		}
	}




	for(int i=0;i<dimx;i++){
		for(int j=0;j<dimy;j++){
			if(partitions[i][j].status==BORDER||partitions[i][j].status==IN){
				cout<<i<<" "<<j;
				if(partitions[i][j].top){
					cout<<" t";
				}
				if(partitions[i][j].bottom){
					cout<<" b";
				}
				if(partitions[i][j].left){
					cout<<" l";
				}
				if(partitions[i][j].right){
					cout<<" r";
				}
				cout<<endl;
				MyPolygon *m = gen_box(start_x+i*step_x,start_y+j*step_y,start_x+(i+1)*step_x,start_y+(j+1)*step_y);
				parts.push_back(m);
			}
		}
	}
	return parts;
}

MyPolygon *MyPolygon::gen_box(double min_x,double min_y,double max_x,double max_y){
	MyPolygon *mbb = new MyPolygon();
	mbb->boundary = new double[11];
	mbb->boundary[0] = 5;
	mbb->boundary[1] = min_x;
	mbb->boundary[2] = min_y;
	mbb->boundary[3] = max_x;
	mbb->boundary[4] = min_y;
	mbb->boundary[5] = max_x;
	mbb->boundary[6] = max_y;
	mbb->boundary[7] = min_x;
	mbb->boundary[8] = max_y;
	mbb->boundary[9] = min_x;
	mbb->boundary[10] = min_y;
	return mbb;
}


MyPolygon *MyPolygon::getMBB(){
	if(this->mbb){
		return mbb;
	}

	double min_x = 180, min_y = 180, max_x = -180, max_y = -180;
	double *cur_b = boundary;
	int num_vertices = (int)*cur_b;
	cur_b++;
	for(int i=0;i<num_vertices;i++){
		if(min_x>*cur_b){
			min_x = *cur_b;
		}
		if(max_x<*cur_b){
			max_x = *cur_b;
		}
		cur_b++;
		if(min_y>*cur_b){
			min_y = *cur_b;
		}
		if(max_y<*cur_b){
			max_y = *cur_b;
		}
		cur_b++;
	}
	mbb = gen_box(min_x,min_y,max_x,max_y);
	return mbb;
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



