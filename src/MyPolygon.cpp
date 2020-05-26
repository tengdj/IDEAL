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
	polygon->getMBB();
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
	polygon->getMBB();
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


vector<vector<Pixel>> MyPolygon::partition(const int dimx, const int dimy){
	assert(dimx>0&&dimy>0);

	if(partitions.size()==dimx){
		if(partitions[0].size()==dimy){
			return partitions;
		}
	}
	for(vector<Pixel> &rows:partitions){
		rows.clear();
	}
	partitions.clear();

	getMBB();
	const double start_x = mbb->low[0];
	const double start_y = mbb->low[1];
	const double end_x = mbb->high[0];
	const double end_y = mbb->high[1];
	step_x = (end_x-start_x)/dimx;
	step_y = (end_y-start_y)/dimy;
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
		int miny = min(cur_starty,cur_endy);
		int maxy = max(cur_starty,cur_endy);

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
				//check horizontally
				if(x!=cur_endx){
					double xval = 0;
					if(cur_startx<cur_endx){
						xval = ((double)x+1)*step_x+start_x;
					}else{
						xval = (double)x*step_x+start_x;
					}
					double yval = xval*a+b;
					int cur_y = (yval-start_y)/step_y;
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
					double yval = 0;
					if(cur_starty<cur_endy){
						yval = (y+1)*step_y+start_y;
					}else{
						yval = y*step_y+start_y;
					}
					double xval = (yval-b)/a;
					int cur_x = (xval-start_x)/step_x;
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
				assert(passed);
			}
		}
	}

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
	return partitions;
}

vector<vector<Pixel>> MyPolygon::decode(char *data){
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

char *MyPolygon::encode(vector<vector<Pixel>> partitions){
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

	int nvert = (int)(*this->boundary);
	bool ret = false;
	for (int i = 0, j = nvert-1; i < nvert; j = i++) {
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




Pixel *MyPolygon::getMBB(){
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
	mbb = new Pixel();
	mbb->low[0] = min_x;
	mbb->low[1] = min_y;
	mbb->high[0] = max_x;
	mbb->high[1] = max_y;
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
	if(false){
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
	}
	cout<<")";
}

void MyPolygon::print(){
	cout<<"POLYGON";
	print_without_head();
	cout<<endl;
}




