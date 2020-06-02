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

vector<MyPolygon *> MyPolygon::load_binary_file(const char *path){
	log("loading polygon from %s",path);
	ifstream infile;
	infile.open(path, ios::in | ios::binary);
	vector<MyPolygon *> polygons;
	int id = 0;
	while(!infile.eof()){
		MyPolygon *poly = read_polygon_binary_file(infile);
		if(!poly){
			continue;
		}
		if(poly->area()<=0){
			delete poly;
			continue;
		}
		poly->setid(id++);
		polygons.push_back(poly);
		if(id%1000000==0){
			log("read %d polygons",id);
		}
	}
	infile.close();
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

void MyPolygon::evaluate_border(int &dimx, int &dimy){
	// normalize
	getMBB();
	const double start_x = mbb->low[0];
	const double start_y = mbb->low[1];
	const double end_x = mbb->high[0];
	const double end_y = mbb->high[1];
	step_x = (end_x-start_x)/dimx;
	step_y = (end_y-start_y)/dimy;

	if(step_x<0.00001){
		step_x = 0.00001;
		dimx = (end_x-start_x)/step_x+1;
	}
	if(step_y<0.00001){
		step_y = 0.00001;
		dimy = (end_y-start_y)/step_y+1;
	}

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

	for(int i=0;i<get_num_vertices()-1;i++){
		double x1 = getx(i);
		double y1 = gety(i);
		double x2 = getx(i+1);
		double y2 = gety(i+1);

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
		if(cur_startx>=dimx||cur_endx>=dimx||cur_starty>=dimy||cur_endy>=dimy){
			cout<<"xrange\t"<<cur_startx<<" "<<cur_endx<<endl;
			cout<<"yrange\t"<<cur_starty<<" "<<cur_endy<<endl;
			printf("xrange_val\t%f %f\n",(x1-start_x)/step_x, (x2-start_x)/step_x);
			printf("yrange_val\t%f %f\n",(y1-start_y)/step_y, (y2-start_y)/step_y);
		}
		assert(cur_startx<dimx);
		assert(cur_endx<dimx);
		assert(cur_starty<dimy);
		assert(cur_endy<dimy);


		partitions[cur_startx][cur_starty].status = BORDER;
		partitions[cur_endx][cur_endy].status = BORDER;

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
			}else { //border[bottom] down
				for(int y=cur_starty;y>cur_endy;y--){
					partitions[cur_startx][y].leave(x1, BOTTOM);
					partitions[cur_startx][y-1].enter(x1, TOP);
				}
			}
		}else{
			// solve the line function
			double a = (y1-y2)/(x1-x2);
			double b = (x1*y2-x2*y1)/(x1-x2);

			int x = cur_startx;
			int y = cur_starty;
			while(x!=cur_endx||y!=cur_endy){
				bool passed = false;
				double yval = 0;
				double xval = 0;
				int cur_x = 0;
				int cur_y = 0;
				//check horizontally
				if(x!=cur_endx){
					if(cur_startx<cur_endx){
						xval = ((double)x+1)*step_x+start_x;
					}else{
						xval = (double)x*step_x+start_x;
					}
					yval = xval*a+b;
					cur_y = (yval-start_y)/step_y;
					//printf("y %f %d\n",(yval-start_y)/step_y,cur_y);
					if(cur_y>max(cur_endy, cur_starty)){
						cur_y=max(cur_endy, cur_starty);
					}
					if(cur_y<min(cur_endy, cur_starty)){
						cur_y=min(cur_endy, cur_starty);
					}
					if(cur_y==y){
						passed = true;
						// left to right
						if(cur_startx<cur_endx){
							partitions[x++][y].leave(yval,RIGHT);
							partitions[x][y].enter(yval,LEFT);
						}else{//right to left
							partitions[x--][y].leave(yval,LEFT);
							partitions[x][y].enter(yval,RIGHT);
						}
					}
				}
				//check vertically
				if(y!=cur_endy){
					if(cur_starty<cur_endy){
						yval = (y+1)*step_y+start_y;
					}else{
						yval = y*step_y+start_y;
					}
					xval = (yval-b)/a;
					int cur_x = (xval-start_x)/step_x;
					//printf("x %f %d\n",(xval-start_x)/step_x,cur_x);
					if(cur_x>max(cur_endx, cur_startx)){
						cur_x=max(cur_endx, cur_startx);
					}
					if(cur_x<min(cur_endx, cur_startx)){
						cur_x=min(cur_endx, cur_startx);
					}
					if(cur_x==x){
						passed = true;
						if(cur_starty<cur_endy){// bottom up
							partitions[x][y++].leave(xval, TOP);
							partitions[x][y].enter(xval, BOTTOM);
						}else{// top down
							partitions[x][y--].leave(xval, BOTTOM);
							partitions[x][y].enter(xval, TOP);
						}
					}
				}
				if(!passed){
					this->print(true);
					cout<<"dim\t"<<dimx<<" "<<dimy<<endl;
					printf("val\t%f %f\n",(xval-start_x)/step_x, (yval-start_y)/step_y);
					cout<<"curxy\t"<<x<<" "<<y<<endl;
					cout<<"calxy\t"<<cur_x<<" "<<cur_y<<endl;
					cout<<"xrange\t"<<cur_startx<<" "<<cur_endx<<endl;
					cout<<"yrange\t"<<cur_starty<<" "<<cur_endy<<endl;
					printf("xrange_val\t%f %f\n",(x1-start_x)/step_x, (x2-start_x)/step_x);
					printf("yrange_val\t%f %f\n",(y1-start_y)/step_y, (y2-start_y)/step_y);
				}
				assert(passed);
			}
		}
	}
}

vector<vector<Pixel>> MyPolygon::partition(int dimx, int dimy){
	assert(dimx>0&&dimy>0);
	pthread_mutex_lock(&partition_lock);
	if(partitioned){
		pthread_mutex_unlock(&partition_lock);
		return partitions;
	}

	evaluate_border(dimx,dimy);

	for(vector<Pixel> &parts:partitions){
		for(Pixel &p:parts){
			p.process_enter_leave();
		}
	}
	//spread the in/out with edge information
	for(int i=0;i<dimx;i++){
		for(int j=0;j<dimy;j++){
			// bottom up
			if(partitions[i][j].border[TOP]==IN&&partitions[i][j+1].status==OUT){
				partitions[i][j+1].setin();
			}else if(partitions[i][j+1].border[BOTTOM]==IN&&partitions[i][j].status==OUT){
				partitions[i][j].setin();
			}
			// left to right
			if(partitions[i][j].border[RIGHT]==IN&&partitions[i+1][j].status==OUT){
				partitions[i+1][j].setin();
			}else if(partitions[i+1][j].border[LEFT]==IN&&partitions[i][j].status==OUT){
				partitions[i][j].setin();
			}
		}
	}
	partitioned = true;
	pthread_mutex_unlock(&partition_lock);
	return partitions;
}



vector<vector<Pixel>> MyPolygon::partition_scanline(int dimx, int dimy){
	assert(dimx>0&&dimy>0);
	pthread_mutex_lock(&partition_lock);
	if(partitioned){
		pthread_mutex_unlock(&partition_lock);
		return partitions;
	}
	//
	evaluate_border(dimx, dimy);
	for(vector<Pixel> &rows:partitions){
		for(Pixel &p:rows){
			if(p.crosses.size()>0){
				p.status = BORDER;
			}
		}
	}
	for(int y=1;y<dimy;y++){
		bool isin = false;
		for(int x=0;x<dimx;x++){
			if(partitions[x][y].status!=BORDER){
				if(isin){
					partitions[x][y].status = IN;
				}
				continue;
			}
			for(cross_info &c:partitions[x][y].crosses){
				if(c.direction==BOTTOM){
					isin = !isin;
				}
			}
		}
	}

	partitioned = true;
	pthread_mutex_unlock(&partition_lock);
	return partitions;
}


vector<vector<Pixel>> MyPolygon::partition_with_query(int dimx, int dimy){
	assert(dimx>0&&dimy>0);
	pthread_mutex_lock(&partition_lock);
	if(partitioned){
		pthread_mutex_unlock(&partition_lock);
		return partitions;
	}
	// normalize
	getMBB();
	const double start_x = mbb->low[0];
	const double start_y = mbb->low[1];
	const double end_x = mbb->high[0];
	const double end_y = mbb->high[1];
	step_x = (end_x-start_x)/dimx;
	step_y = (end_y-start_y)/dimy;

	if(step_x<0.00001){
		step_x = 0.00001;
		dimx = (end_x-start_x)/step_x+1;
	}
	if(step_y<0.00001){
		step_y = 0.00001;
		dimy = (end_y-start_y)/step_y+1;
	}

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
			Point p1(m.low[0],m.low[1]);
			Point p2(m.high[0],m.high[1]);
			bool in1 = contain(p1,false);
			bool in2 = contain(p2,false);;
			if(in1&&in2){
				m.status = IN;
			}else if(!in1&&!in2){
				m.status = OUT;
			}else{
				m.status = BORDER;
			}
		}
		partitions.push_back(v);
	}

	partitioned = true;
	pthread_mutex_unlock(&partition_lock);
	return partitions;
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

			p.status = (PartitionStatus)((cur>>shift)&3);
			p.id[0] = i;
			p.id[1] = j;
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

/*
 * query
 *
 * */

bool MyPolygon::contain(Point &p,bool use_partition){

	if(!mbb){
		getMBB();
	}
	if(mbb->low[0]>p.x||mbb->low[1]>p.y||
	   mbb->high[0]<p.x||mbb->high[1]<p.y){
		return false;
	}

	if(partitioned&&use_partition){
		int part_x = this->get_pixel_x(p.x);
		int part_y = this->get_pixel_y(p.y);
		if(partitions[part_x][part_y].status==IN){
			return true;
		}
		if(partitions[part_x][part_y].status==OUT){
			return false;
		}
	}

	bool ret = false;
	for (int i = 0, j = get_num_vertices()-1; i < get_num_vertices(); j = i++) {
		// segment i->j intersect with line y=p.y
		if ( ((gety(i)>p.y) != (gety(j)>p.y))){
			double a = (getx(j)-getx(i)) / (gety(j)-gety(i));
			if(p.x - getx(i) < a * (p.y-gety(i))){
				ret = !ret;
			}
		}
	}
	return ret;

}


bool MyPolygon::contain_try_partition(MyPolygon *target, query_context *ctx){
	assert(partitioned&&ctx->use_partition);
	ctx->partition_determined = true;
	Pixel *a = getMBB();
	Pixel *b = target->getMBB();
	if(!a->contain(b)){
		return false;
	}
	// test all the pixels
	int txstart = this->get_pixel_x(a->low[0]);
	int txend = this->get_pixel_x(a->high[0]);
	int tystart = this->get_pixel_y(a->low[1]);
	int tyend = this->get_pixel_y(a->high[0]);
	int incount = 0;
	for(int i=txstart;i<=txend;i++){
		for(int j=tystart;j<=tyend;j++){
			if(partitions[i][j].status==OUT){
				return false;
			}
			if(partitions[i][j].status==IN){
				incount++;
			}
		}
	}
	// all is in
	if(incount==(txend-txstart+1)*(tyend-tystart+1)){
		return true;
	}

	ctx->partition_determined = false;
	return false;
}

bool MyPolygon::contain(MyPolygon *target, query_context *ctx){
	// check each vertex in target
	for(int i=0;i<target->get_num_vertices();i++){
		Point p(target->getx(i),target->gety(i));
		if(!contain(p, ctx->use_partition)){
			return false;
		}
	}
	return true;
}

bool MyPolygon::contain(Pixel *target, bool use_partition){

	Pixel *a = getMBB();
	if(!a->contain(target)){
		return false;
	}
	if(partitioned&&use_partition){
		// test all the pixels
		int txstart = this->get_pixel_x(a->low[0]);
		int txend = this->get_pixel_x(a->high[0]);
		int tystart = this->get_pixel_y(a->low[1]);
		int tyend = this->get_pixel_y(a->high[0]);
		int incount = 0;
		for(int i=txstart;i<=txend;i++){
			for(int j=tystart;j<=tyend;j++){
				if(partitions[i][j].status==OUT){
					return false;
				}
				if(partitions[i][j].status==IN){
					incount++;
				}
			}
		}
		// all is in
		if(incount==(txend-txstart+1)*(tyend-tystart+1)){
			return true;
		}
	}

	Point p1(target->low[0], target->low[1]);
	Point p2(target->high[0], target->high[1]);
	return contain(p1, use_partition)&&contain(p2, use_partition);

}

bool MyPolygon::intersect(MyPolygon *target, bool use_partition){

	Pixel *a = getMBB();
	Pixel *b = target->getMBB();
	if(!a->intersect(b)){
		return false;
	}
	if(partitioned&&use_partition){
		// test all the pixels
		int txstart = this->get_pixel_x(a->low[0]);
		int txend = this->get_pixel_x(a->high[0]);
		int tystart = this->get_pixel_y(a->low[1]);
		int tyend = this->get_pixel_y(a->high[0]);
		int outcount = 0;
		for(int i=txstart;i<=txend;i++){
			for(int j=tystart;j<=tyend;j++){
				if(partitions[i][j].status==OUT){
					outcount++;
				}else if(partitions[i][j].status==IN){
					return true;
				}
			}
		}
		// all is out
		if(outcount==(txend-txstart+1)*(tyend-tystart+1)){
			return false;
		}
	}

	// check each vertex in target
	for(int i=0;i<target->get_num_vertices();i++){
		Point p(target->getx(i),target->gety(i));
		if(contain(p, use_partition)){
			return true;
		}
	}
	return false;
}




