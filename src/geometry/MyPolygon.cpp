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
#include <thread>

using namespace std;

VertexSequence::VertexSequence(int nv){
	x = new double[nv];
	y = new double[nv];
	num_vertices = nv;
}
VertexSequence::VertexSequence(int nv, double *xx, double *yy){
	num_vertices = nv;
	x = new double[nv];
	y = new double[nv];
	memcpy((char *)x,(char *)xx,num_vertices*sizeof(double));
	memcpy((char *)y,(char *)yy,num_vertices*sizeof(double));
};

VertexSequence::~VertexSequence(){
	if(x){
		delete []x;
	}
	if(y){
		delete []y;
	}
}

vector<Point *> VertexSequence::pack_to_polyline(){
	vector<Point *> polyline;
	polyline.resize(num_vertices-1);
	for(int i=0;i<num_vertices-1;i++){
		polyline[i] = new Point(x[i],y[i]);
	}
	return polyline;
}

Pixel *VertexSequence::getMBR(){
	Pixel *mbr = new Pixel();
	double min_x = 180, min_y = 180, max_x = -180, max_y = -180;
	for(int i=0;i<num_vertices;i++){
		if(min_x>x[i]){
			min_x = x[i];
		}
		if(max_x<x[i]){
			max_x = x[i];
		}
		if(min_y>y[i]){
			min_y = y[i];
		}
		if(max_y<y[i]){
			max_y = y[i];
		}
	}
	mbr->low[0] = min_x;
	mbr->low[1] = min_y;
	mbr->high[0] = max_x;
	mbr->high[1] = max_y;
	return mbr;
}

void VertexSequence::fix(){
	if(num_vertices<=3){
		return;
	}
	int cur = 0;
	int cur_store = 0;
	while(cur<num_vertices-1){
		if(!double_equal(x[cur],x[cur+1])||!double_equal(y[cur],y[cur+1])){
			x[cur_store] = x[cur];
			y[cur_store] = y[cur];
			cur_store++;
		}
		cur++;
	}
	x[cur_store] = x[cur];
	y[cur_store] = y[cur];
	num_vertices = cur_store+1;
}

size_t VertexSequence::encode(char *dest){
	size_t encoded = 0;
	((long *)dest)[0] = num_vertices;
	encoded += sizeof(long);
	memcpy(dest+encoded,(char *)x,num_vertices*sizeof(double));
	encoded += num_vertices*sizeof(double);
	memcpy(dest+encoded,(char *)y,num_vertices*sizeof(double));
	encoded += num_vertices*sizeof(double);
	return encoded;
}

size_t VertexSequence::decode(char *source){
	size_t decoded = 0;
	num_vertices = ((long *)source)[0];
	x = new double[num_vertices];
	y = new double[num_vertices];
	decoded += sizeof(long);
	memcpy((char *)x,source+decoded,num_vertices*sizeof(double));
	decoded += num_vertices*sizeof(double);
	memcpy((char *)y,source+decoded,num_vertices*sizeof(double));
	decoded += num_vertices*sizeof(double);
	return decoded;
}

size_t VertexSequence::get_data_size(){
	return sizeof(long)+num_vertices*2*sizeof(double);
}


VertexSequence *VertexSequence::clone(){
	VertexSequence *ret = new VertexSequence(num_vertices);
	memcpy((void *)ret->x,(void *)x,sizeof(double)*num_vertices);
	memcpy((void *)ret->y,(void *)y,sizeof(double)*num_vertices);
	return ret;
}
void VertexSequence::print(){
	cout<<"(";
	for(int i=0;i<num_vertices;i++){
		if(i!=0){
			cout<<",";
		}
		printf("%f ",x[i]);
		printf("%f",y[i]);
	}
	cout<<")";
}

double VertexSequence::area(){
	double sum = 0;
	for(int i=0;i<num_vertices-1;i++){
		sum += (x[i]-x[i+1])*(y[i+1]+y[i]);
	}
	return sum;
}

bool VertexSequence::clockwise(){
	return area()<0;
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
	VertexSequence *vs = new VertexSequence(num_vertices);

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
	if(poly->boundary->clockwise()){
		poly->boundary->reverse();
	}
//	string str = poly->to_string();
	poly->boundary->fix();
//	if(poly->boundary->num_vertices<10){
//		cout<<str<<endl;
//		poly->print();
//	}
	for(int i=0;i<num_holes;i++){
		infile.read((char *)&num_vertices,sizeof(long));
		assert(num_vertices);
		VertexSequence *vs = new VertexSequence(num_vertices);
		infile.read((char *)vs->x,num_vertices*sizeof(double));
		infile.read((char *)vs->y,num_vertices*sizeof(double));
		if(!vs->clockwise()){
			vs->reverse();
		}
		vs->fix();
		poly->internal_polygons.push_back(vs);
	}
	return poly;
}

inline bool compareIterator(MyPolygon *p1, MyPolygon *p2){
	return p1->get_num_vertices()>p2->get_num_vertices();
}

MyPolygon * MyPolygon::load_binary_file_single(const char *path, query_context ctx, int idx){
	ifstream infile;
	infile.open(path, ios::in | ios::binary);
	size_t num;
	infile.read((char *)&num, sizeof(size_t));
	size_t offset = sizeof(size_t)+idx*sizeof(size_t);
	infile.seekg(offset,infile.beg);
	infile.read((char *)&offset, sizeof(size_t));
	//cout<<num<<" "<<offset<<endl;

	infile.seekg(offset,infile.beg);
	MyPolygon *poly = read_polygon_binary_file(infile);
	infile.close();
	assert(poly);
	return poly;

}


vector<MyPolygon *> MyPolygon::load_binary_file(const char *path, query_context ctx, bool sample){
	vector<MyPolygon *> polygons;
	if(!file_exist(path)){
		log("%s does not exist",path);
		return polygons;
	}
	log("loading polygon from %s",path);
	struct timeval start = get_cur_time();
	ifstream infile;
	infile.open(path, ios::in | ios::binary);
	size_t off;
	size_t num_polygons = 0;
	//seek to the first polygon
	infile.read((char *)&num_polygons, sizeof(size_t));
	size_t *offsets = new size_t[num_polygons];
	infile.read((char *)offsets, sizeof(size_t)*num_polygons);
	size_t num_edges = 0;
	size_t data_size = 0;
	int next = 10;
	for(int i=0;i<num_polygons;i++){
		infile.seekg(offsets[i], infile.beg);
		MyPolygon *poly = read_polygon_binary_file(infile);
		assert(poly);
		num_edges += poly->get_num_vertices();
		data_size += poly->get_data_size();
		poly->setid(i);
		polygons.push_back(poly);
		if(i*100/num_polygons>=next){
			log("loaded %d%%",next);
			next+=10;
		}
	}
	infile.close();
	if(ctx.sort_polygons){
		std::sort(polygons.begin(),polygons.end(),compareIterator);
	}

	logt("loaded %ld polygons each with %ld edges %ld MB", start, polygons.size(),num_edges/polygons.size(),data_size/1024/1024);
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
	polygon->boundary = read_vertices(wkt, offset,false);
	if(polygon->boundary->clockwise()){
		polygon->boundary->reverse();
	}
	polygon->boundary->fix();
	skip_space(wkt, offset);
	//polygons as the holes of the boundary polygon
	while(wkt[offset]==','){
		offset++;
		VertexSequence *vc = read_vertices(wkt, offset,true);
		if(!vc->clockwise()){
			vc->reverse();
		}
		vc->fix();
		polygon->internal_polygons.push_back(vc);

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
	if(mbr){
		delete mbr;
	}

	if(mer){
		delete mer;
	}
	if(raster){
		delete raster;
	}
	if(qtree){
		delete qtree;
	}
	if(rtree){
		delete rtree;
	}
	for(Point *p:polyline){
		delete p;
	}
	polyline.clear();
	for(Triangle *tri:triangles){
		delete tri;
	}
	triangles.clear();
}


size_t MyPolygon::get_data_size(){
	size_t ds = 0;
	ds += sizeof(long);
	// for boundary
	ds += boundary->get_data_size();
//	for(VertexSequence *vs:internal_polygons){
//		ds += vs->get_data_size();
//	}
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
	MyPolygon *mbr = new MyPolygon();
	mbr->boundary = new VertexSequence(5);
	mbr->boundary->x[0] = min_x;
	mbr->boundary->y[0] = min_y;
	mbr->boundary->x[1] = max_x;
	mbr->boundary->y[1] = min_y;
	mbr->boundary->x[2] = max_x;
	mbr->boundary->y[2] = max_y;
	mbr->boundary->x[3] = min_x;
	mbr->boundary->y[3] = max_y;
	mbr->boundary->x[4] = min_x;
	mbr->boundary->y[4] = min_y;
	return mbr;
}

MyPolygon *MyPolygon::gen_box(Pixel &pix){
	return gen_box(pix.low[0],pix.low[1],pix.high[0],pix.high[1]);
}

vector<Point> MyPolygon::generate_test_points(int num){
	getMBB();
	vector<Point> points;
	for(int i=0;i<num;i++){
		Point p;
		p.x = get_rand_double()*(mbr->high[0]-mbr->low[0])+mbr->low[0];
		p.y = get_rand_double()*(mbr->high[1]-mbr->low[1])+mbr->low[1];
		points.push_back(p);
	}
	return points;
}


vector<MyPolygon *> MyPolygon::generate_test_polygons(int num){
	getMBB();
	vector<MyPolygon *> polys;
	for(int i=0;i<num;i++){
		double start_x = get_rand_double()*(mbr->high[0]-mbr->low[0])+mbr->low[0];
		double start_y = get_rand_double()*(mbr->high[1]-mbr->low[1])+mbr->low[1];

		double end_x = get_rand_double()*(mbr->high[0]-mbr->low[0])/100+start_x;
		double end_y = get_rand_double()*(mbr->high[1]-mbr->low[1])/100+start_y;
		MyPolygon *m = gen_box(start_x, start_y, end_x, end_y);
		polys.push_back(m);
	}

	return polys;
}


void MyPolygon::print_without_head(bool print_hole){
	assert(boundary);
	cout<<"(";

	boundary->print();
	if(print_hole){
		for(VertexSequence *vs:this->internal_polygons){
			cout<<", ";
			vs->print();
		}
	}
	cout<<")";
}

void MyPolygon::print_triangles(){
	printf("MULTIPOLYGON(");
	bool first = true;
	for(Triangle *tri:triangles){
		if(!first){
			printf(",");
		}else{
			first = false;
		}
		tri->print();
	}
	printf(")\n");
}

void MyPolygon::print(bool print_hole){
	cout<<"id:\t"<<this->id<<endl;
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
				MyPolygon *m = gen_box(cur->mbr);
				if(cur->interior){
					inpolys->insert_polygon(m);
				}else if(cur->exterior){
					outpolys->insert_polygon(m);
				}else{
					borderpolys->insert_polygon(m);
				}
			}else{
				cur->push(ws);
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

void MyPolygon::rasterization(int vpr){
	assert(vpr>0);

	pthread_mutex_lock(&ideal_partition_lock);
	if(raster==NULL){
		raster = new MyRaster(boundary,vpr);
		raster->rasterization();
	}
	pthread_mutex_unlock(&ideal_partition_lock);
}


QTNode *MyPolygon::partition_qtree(const int vpr){

	pthread_mutex_lock(&qtree_partition_lock);
	if(qtree){
		pthread_mutex_unlock(&qtree_partition_lock);
		return qtree;
	}

	int num_boxes = get_num_vertices()/vpr;
	int level = 1;
	while(pow(4,level)<num_boxes){
		level++;
	}
	int dimx = pow(2,level);
	int dimy = dimx;

	MyRaster *ras = new MyRaster(boundary,dimx,dimy);
	ras->rasterization();
	int box_count = 4;
	int cur_level = 1;

	qtree = new QTNode(*(this->getMBB()));
	std::stack<QTNode *> ws;
	qtree->split();
	qtree->push(ws);

	vector<QTNode *> level_nodes;

	//breadth first traverse
	query_context qc;
	while(box_count<num_boxes||!ws.empty()){
		if(cur_level>level){
			dimx *= 2;
			dimy *= 2;
			level = cur_level;
			delete ras;
			ras = new MyRaster(boundary,dimx,dimy);
			ras->rasterization();
		}

		while(!ws.empty()){
			QTNode *cur = ws.top();
			ws.pop();
			bool contained = false;
			if(!ras->contain(&cur->mbr, contained)){
				level_nodes.push_back(cur);
			}else if(contained){
				cur->interior = true;
			}else{
				cur->exterior = true;
			}
		}


		for(QTNode *n:level_nodes){
			if(box_count<num_boxes){
				n->split();
				n->push(ws);
				box_count += 3;
			}else{
				break;
			}
		}
		cur_level++;
		level_nodes.clear();
	}
	//assert(box_count>=num_boxes);

	delete ras;
	pthread_mutex_unlock(&qtree_partition_lock);

	return qtree;
}



vector<vector<Pixel>> MyPolygon::decode_partition(char *data){
	assert(data);
	vector<vector<Pixel>> partitions;

	int dimx = data[0];
	int dimy = data[1];
	if(dimx==0||dimy==0||dimx>255||dimy>255){
		return partitions;
	}
	double *meta = (double *)(data+2);
	double start_x = meta[0], start_y = meta[1], step_x = meta[2], step_y = meta[3];

	int index = 2+4*8;
	for(int i=0;i<dimx;i++){
		char cur = data[index++];
		int shift = 0;
		vector<Pixel> row;
		for(int j=0;j<dimy;j++){
			Pixel p;
			p.id[0] = i;
			p.id[1] = j;
			p.status = (PartitionStatus)((cur>>shift)&3);
			p.low[0] = i*step_x+start_x;
			p.low[1] = j*step_y+start_y;
			p.high[0] = (i+1)*step_x+start_x;
			p.high[1] = (j+1)*step_y+start_y;
			row.push_back(p);
			if(shift==6){
				cur = data[index++];
				shift=0;
			}else{
				shift+=2;
			}
		}
		partitions.push_back(row);
	}

	return partitions;
}

char *MyPolygon::encode_partition(vector<vector<Pixel>> partitions){
	int dimx = partitions.size();
	if(dimx==0){
		return NULL;
	}
	int dimy = partitions[0].size();
	if(dimy==0){
		return NULL;
	}
	for(vector<Pixel> &rows:partitions){
		 if(rows.size()!=dimy){
			 return NULL;
		 }
	}
	assert(dimx<=255&&dimy<=255);

	char *data = new char[((dimy+3)/4*dimx)+2+4*8];
	data[0] = (char)dimx;
	data[1] = (char)dimy;
	double *meta = (double *)(data+2);
	meta[0] = partitions[0][0].low[0];
	meta[1] = partitions[0][0].low[1];
	meta[2] = partitions[0][0].high[0]-partitions[0][0].low[0];
	meta[3] = partitions[0][0].high[1]-partitions[0][0].low[1];
	int index = 2+4*8;
	for(int i=0;i<dimx;i++){
		char cur = 0;
		int shift = 0;
		for(int j=0;j<dimy;j++){
			cur |= ((uint8_t)partitions[i][j].status)<<shift;
			if(shift==6){
				data[index++]=cur;
				cur = 0;
				shift=0;
			}else{
				shift+=2;
			}
		}
		if(shift>0){
			data[index++]=cur;
		}
	}
	return data;
}



size_t MyPolygon::partition_size(){
	size_t size = 0;
	const int nump = raster->get_num_pixels();

	const int nump_border = raster->get_num_pixels(BORDER);
	int bits_b = 0;
	int tmp = nump_border+2;
	while(tmp>0){
		bits_b++;
		tmp /= 2;
	}

	int numv = this->get_num_vertices();
	tmp = numv;
	int bits_v = 0;
	while(tmp>0){
		bits_v++;
		tmp /= 2;
	}

	bits_v = bits_v*2;

//	bits_b = (2+bits_b+7)/8*8;
//	bits_v = (2*bits_v+7)/8*8;

	int numc = raster->get_num_crosses();
	return (bits_b*nump + bits_v*nump_border + numc*64+7)/8;
}






