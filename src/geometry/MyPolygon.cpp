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
#include <fstream>
#include <float.h>
#include <sstream>
#include <vector>

using namespace std;

double VertexSequence::area(){
	double sum = 0;
	for(int i=0;i<num_vertices-1;i++){
		sum += (x[i+1]-x[i])*(y[i+1]+y[i]);
	}
	return sum;
}

bool VertexSequence::clockwise(){
	return area()>0;
}

void VertexSequence::reverse(){
	double tmp = 0;
	for(int i=0;i<num_vertices/2;i++){
		tmp = x[i];
		x[i] = x[num_vertices-1-i];
		x[num_vertices-1-i] = tmp;
		tmp = y[i];
		y[i] = y[num_vertices-1-i];
		y[num_vertices-1-i] = tmp;
	}
}

VertexSequence *MyPolygon::read_vertices(const char *wkt, size_t &offset, bool clockwise){
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
	VertexSequence *vs = new VertexSequence();
	vs->num_vertices = num_vertices;
	vs->x = new double[num_vertices];
	vs->y = new double[num_vertices];

	// read x/y
	for(int i=0;i<num_vertices;i++){
		vs->x[i] = read_double(wkt, offset);
		vs->y[i] = read_double(wkt, offset);
	}
	if(clockwise){
		if(!vs->clockwise()){
			vs->reverse();
		}
	}else{
		if(vs->clockwise()){
			vs->reverse();
		}
	}

	// read until the right parenthesis
	skip_space(wkt, offset);
	assert(wkt[offset++]==')');
	return vs;
}


MyPolygon * MyPolygon::read_polygon_binary_file(ifstream &infile){
	long num_holes = 0;
	infile.read((char *)&num_holes,sizeof(long));
	long num_vertices = 0;
	infile.read((char *)&num_vertices,sizeof(long));
	if(num_vertices==0){
		return NULL;
	}
	MyPolygon *poly = new MyPolygon();
	poly->boundary = new VertexSequence(num_vertices);
	infile.read((char *)poly->boundary->x,num_vertices*sizeof(double));
	infile.read((char *)poly->boundary->y,num_vertices*sizeof(double));
	for(int i=0;i<num_holes;i++){
		infile.read((char *)&num_vertices,sizeof(long));
		assert(num_vertices);
		VertexSequence *vs = new VertexSequence(num_vertices);
		infile.read((char *)vs->x,num_vertices*sizeof(double));
		infile.read((char *)vs->y,num_vertices*sizeof(double));
		poly->internal_polygons.push_back(vs);
	}
	return poly;
}

inline bool compareIterator(MyPolygon *p1, MyPolygon *p2){
	return p1->get_num_vertices()>p2->get_num_vertices();
}

vector<MyPolygon *> MyPolygon::load_binary_file(const char *path, query_context ctx){
	vector<MyPolygon *> polygons;
	if(!file_exist(path)){
		log("%s does not exist",path);
		return polygons;
	}
	log("loading polygon from %s",path);
	ifstream infile;
	infile.open(path, ios::in | ios::binary);
	int id = 0;
	int iii = 0;
	while(!infile.eof()){
		MyPolygon *poly = read_polygon_binary_file(infile);
		if(!poly){
			continue;
		}

		if(poly->get_num_vertices()<ctx.small_threshold||
			  poly->get_num_vertices()>ctx.big_threshold||
				poly->area()<=0){
			delete poly;
			continue;
		}
		if(ctx.sample_rate<1.0 && !tryluck(ctx.sample_rate)){
			if(poly){
				delete poly;
			}
			continue;
		}
		poly->setid(id++);
		polygons.push_back(poly);
		if(id%1000000==0){
			log("read %d polygons",id);
		}
	}
	infile.close();
	if(ctx.sort_polygons){
		std::sort(polygons.begin(),polygons.end(),compareIterator);
	}
	return polygons;
}


MyPolygon *MyPolygon::clone(){
	MyPolygon *polygon = new MyPolygon();
	polygon->boundary = boundary->clone();
	for(VertexSequence *t:internal_polygons){
		polygon->internal_polygons.push_back(t->clone());
	}
	polygon->getMBB();
	return polygon;
}


MyPolygon *MyPolygon::read_polygon(const char *wkt, size_t &offset){
	MyPolygon *polygon = new MyPolygon();
	skip_space(wkt, offset);
	// left parentheses for the entire polygon
	assert(wkt[offset++]=='(');

	// read the vertices of the boundary polygon
	// the vertex must rotation in clockwise
	polygon->boundary = read_vertices(wkt, offset,true);

	skip_space(wkt, offset);
	//polygons as the holes of the boundary polygon
	while(wkt[offset]==','){
		offset++;
		polygon->internal_polygons.push_back(read_vertices(wkt, offset,false));

		skip_space(wkt, offset);
	}
	assert(wkt[offset++]==')');
	polygon->getMBB();
	return polygon;
}

MyPolygon::~MyPolygon(){
	if(this->boundary){
		delete boundary;
	}
	for(VertexSequence *p:internal_polygons){
		if(p){
			delete p;
		}
	}
	internal_polygons.clear();
	if(mbb){
		delete mbb;
	}
}


size_t MyPolygon::get_data_size(){
	size_t ds = 0;
	ds += sizeof(long);
	// for boundary
	ds += boundary->get_data_size();
	for(VertexSequence *vs:internal_polygons){
		ds += vs->get_data_size();
	}
	return ds;
}

size_t MyPolygon::encode_to(char *target){
	size_t encoded = 0;
	((long *)target)[0] = internal_polygons.size();
	encoded += sizeof(long);
	encoded += boundary->encode(target+encoded);
	for(VertexSequence *vs:internal_polygons){
		encoded += vs->encode(target+encoded);
	}
	return encoded;
}
size_t MyPolygon::decode_from(char *source){
	size_t decoded = 0;
	assert(!boundary);
	boundary = new VertexSequence();
	size_t num_holes = ((long *)source)[0];
	decoded += sizeof(long);
	decoded += boundary->decode(source+decoded);
	for(size_t i=0;i<num_holes;i++){
		VertexSequence *vs = new VertexSequence();
		decoded += vs->decode(source+decoded);
	}
	return decoded;
}

MyPolygon *MyPolygon::gen_box(double min_x,double min_y,double max_x,double max_y){
	MyPolygon *mbb = new MyPolygon();
	mbb->boundary = new VertexSequence(5);
	mbb->boundary->x[0] = min_x;
	mbb->boundary->y[0] = min_y;
	mbb->boundary->x[1] = min_x;
	mbb->boundary->y[1] = max_y;
	mbb->boundary->x[2] = max_x;
	mbb->boundary->y[2] = max_y;
	mbb->boundary->x[3] = max_x;
	mbb->boundary->y[3] = min_y;
	mbb->boundary->x[4] = min_x;
	mbb->boundary->y[4] = min_y;
	return mbb;
}

vector<Point> MyPolygon::generate_test_points(int num){
	getMBB();
	vector<Point> points;
	for(int i=0;i<num;i++){
		Point p;
		p.x = get_rand_double()*(mbb->high[0]-mbb->low[0])+mbb->low[0];
		p.y = get_rand_double()*(mbb->high[1]-mbb->low[1])+mbb->low[1];
		points.push_back(p);
	}
	return points;
}


vector<MyPolygon *> MyPolygon::generate_test_polygons(int num){
	getMBB();
	vector<MyPolygon *> polys;
	for(int i=0;i<num;i++){
		double start_x = get_rand_double()*(mbb->high[0]-mbb->low[0])+mbb->low[0];
		double start_y = get_rand_double()*(mbb->high[1]-mbb->low[1])+mbb->low[1];

		double end_x = get_rand_double()*(mbb->high[0]-mbb->low[0])/100+start_x;
		double end_y = get_rand_double()*(mbb->high[1]-mbb->low[1])/100+start_y;
		MyPolygon *m = gen_box(start_x, start_y, end_x, end_y);
		polys.push_back(m);
	}

	return polys;
}



Pixel *MyPolygon::getMBB(){
	if(this->mbb){
		return mbb;
	}

	double min_x = 180, min_y = 180, max_x = -180, max_y = -180;
	for(int i=0;i<boundary->num_vertices;i++){
		if(min_x>boundary->x[i]){
			min_x = boundary->x[i];
		}
		if(max_x<boundary->x[i]){
			max_x = boundary->x[i];
		}
		if(min_y>boundary->y[i]){
			min_y = boundary->y[i];
		}
		if(max_y<boundary->y[i]){
			max_y = boundary->y[i];
		}
	}
	mbb = new Pixel();
	mbb->low[0] = min_x;
	mbb->low[1] = min_y;
	mbb->high[0] = max_x;
	mbb->high[1] = max_y;
	return mbb;
}

void MyPolygon::print_without_head(bool print_hole){
	assert(boundary);
	cout<<"((";
	for(int i=0;i<boundary->num_vertices;i++){
		if(i!=0){
			cout<<",";
		}
		printf("%f ",boundary->x[i]);
		printf("%f",boundary->y[i]);
	}
	cout<<")";
	if(print_hole){
		for(VertexSequence *vs:this->internal_polygons){
			cout<<", (";
			for(int i=0;i<vs->num_vertices;i++){
				if(i!=0){
					cout<<",";
				}
				printf("%f ",vs->x[i]);
				printf("%f",vs->y[i]);
			}
			cout<<")";
		}
	}
	cout<<")";
}

void MyPolygon::print(bool print_hole){
	cout<<"POLYGON";
	print_without_head(print_hole);
	cout<<endl;
}
void MyPolygon::print_without_return(bool print_hole){
	cout<<"POLYGON";
	print_without_head(print_hole);
}

string MyPolygon::to_string(bool clockwise){
	std::stringstream ss;
	char double_str[200];
	ss<<"POLYGON";
	ss<<"((";
	if(clockwise){
		for(int i=0;i<boundary->num_vertices;i++){
			if(i!=0){
				ss<<",";
			}
			sprintf(double_str,"%f %f",boundary->x[i],boundary->y[i]);
			ss<<double_str;
		}
	}else{
		for(int i=boundary->num_vertices-1;i>=0;i--){
			if(i!=boundary->num_vertices-1){
				ss<<",";
			}
			sprintf(double_str,"%f %f",boundary->x[i],boundary->y[i]);
			ss<<double_str;
		}
	}

	ss<<"))";
	return ss.str();
}


void MyPolygon::print_partition(query_context qt){
	MyMultiPolygon *inpolys = new MyMultiPolygon();
	MyMultiPolygon *borderpolys = new MyMultiPolygon();
	MyMultiPolygon *outpolys = new MyMultiPolygon();

	if(qt.use_qtree){
		partition_qtree(qt.vpr);
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
	}
	if(qt.use_grid){
		partition(qt.vpr);
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

	cout<<"border:"<<endl;
	borderpolys->print();
	cout<<"in:"<<endl;
	inpolys->print();
	cout<<"out:"<<endl;
	outpolys->print();


	delete borderpolys;
	delete inpolys;
	delete outpolys;
}



Point *Point::read_one_point(string &input_line){

	if(input_line.size()==0){
		return NULL;
	}
	const char *wkt = input_line.c_str();
	size_t offset = 0;
	// read the symbol MULTIPOLYGON
	while(wkt[offset]!='P'){
		offset++;
	}
	for(int i=0;i<strlen(point_char);i++){
		assert(wkt[offset++]==point_char[i]);
	}
	skip_space(wkt,offset);
	Point *p = new Point();
	p->x = read_double(wkt,offset);
	p->y = read_double(wkt,offset);

	return p;
}


