/*
 * MyRaster.cpp
 *
 *  Created on: Jan 4, 2021
 *      Author: teng
 */


#include "MyPolygon.h"

/*
 *
 * functions for rasterization
 *
 *
 * */

MyRaster::MyRaster(VertexSequence *vst, int epp){
	assert(epp>0);
	vs = vst;
	mbr = vs->getMBR();

	double multi = abs((mbr->high[1]-mbr->low[1])/(mbr->high[0]-mbr->low[0]));
	dimx = std::pow((vs->num_vertices/epp)/multi,0.5);
	dimy = dimx*multi;

	if(dimx==0){
		dimx = 1;
	}
	if(dimy==0){
		dimy = 1;
	}

	step_x = (mbr->high[0]-mbr->low[0])/dimx;
	step_y = (mbr->high[1]-mbr->low[1])/dimy;

	if(step_x<0.00001){
		step_x = 0.00001;
		dimx = (mbr->high[0]-mbr->low[0])/step_x+1;
	}
	if(step_y<0.00001){
		step_y = 0.00001;
		dimy = (mbr->high[1]-mbr->low[1])/step_y+1;
	}

}


MyRaster::MyRaster(VertexSequence *vst, int dx, int dy){

	vs = vst;

	mbr = vs->getMBR();
	dimx = dx;
	dimy = dy;

	step_x = (mbr->high[0]-mbr->low[0])/dimx;
	step_y = (mbr->high[1]-mbr->low[1])/dimy;

	if(step_x<0.00001){
		step_x = 0.00001;
		dimx = (mbr->high[0]-mbr->low[0])/step_x+1;
	}
	if(step_y<0.00001){
		step_y = 0.00001;
		dimy = (mbr->high[1]-mbr->low[1])/step_y+1;
	}
}

MyRaster::~MyRaster(){
	for(vector<Pixel *> &rows:pixels){
		for(Pixel *p:rows){
			delete p;
		}
		rows.clear();
	}
	pixels.clear();
}

void MyRaster::init_pixels(){
	assert(mbr);
	const double start_x = mbr->low[0];
	const double start_y = mbr->low[1];
	for(double i=0;i<=dimx;i++){
		vector<Pixel *> v;
		for(double j=0;j<=dimy;j++){
			Pixel *m = new Pixel();
			m->id[0] = i;
			m->id[1] = j;
			m->low[0] = i*step_x+start_x;
			m->high[0] = (i+1.0)*step_x+start_x;
			m->low[1] = j*step_y+start_y;
			m->high[1] = (j+1.0)*step_y+start_y;
			v.push_back(m);
		}
		pixels.push_back(v);
	};
}

void MyRaster::evaluate_edges(){
	// normalize
	assert(mbr);
	const double start_x = mbr->low[0];
	const double start_y = mbr->low[1];

	for(int i=0;i<vs->num_vertices-1;i++){
		double x1 = vs->x[i];
		double y1 = vs->y[i];
		double x2 = vs->x[i+1];
		double y2 = vs->y[i+1];

		int cur_startx = (x1-start_x)/step_x;
		int cur_endx = (x2-start_x)/step_x;
		int cur_starty = (y1-start_y)/step_y;
		int cur_endy = (y2-start_y)/step_y;

		if(cur_startx==dimx+1){
			cur_startx--;
		}
		if(cur_endx==dimx+1){
			cur_endx--;
		}

		int minx = min(cur_startx,cur_endx);
		int maxx = max(cur_startx,cur_endx);

		if(cur_starty==dimy+1){
			cur_starty--;
		}
		if(cur_endy==dimy+1){
			cur_endy--;
		}
		// should not happen for normal cases
		if(cur_startx>dimx||cur_endx>dimx||cur_starty>dimy||cur_endy>dimy){
			cout<<"xrange\t"<<cur_startx<<" "<<cur_endx<<endl;
			cout<<"yrange\t"<<cur_starty<<" "<<cur_endy<<endl;
			printf("xrange_val\t%f %f\n",(x1-start_x)/step_x, (x2-start_x)/step_x);
			printf("yrange_val\t%f %f\n",(y1-start_y)/step_y, (y2-start_y)/step_y);
		}
		assert(cur_startx<=dimx);
		assert(cur_endx<=dimx);
		assert(cur_starty<=dimy);
		assert(cur_endy<=dimy);

		pixels[cur_startx][cur_starty]->status = BORDER;
		pixels[cur_endx][cur_endy]->status = BORDER;

		//in the same pixel
		if(cur_startx==cur_endx&&cur_starty==cur_endy){
			continue;
		}

		if(y1==y2){
			//left to right
			if(cur_startx<cur_endx){
				for(int x=cur_startx;x<cur_endx;x++){
					pixels[x][cur_starty]->leave(y1,RIGHT,i);
					pixels[x+1][cur_starty]->enter(y1,LEFT,i);
				}
			}else { // right to left
				for(int x=cur_startx;x>cur_endx;x--){
					pixels[x][cur_starty]->leave(y1, LEFT,i);
					pixels[x-1][cur_starty]->enter(y1, RIGHT,i);
				}
			}
		}else if(x1==x2){
			//bottom up
			if(cur_starty<cur_endy){
				for(int y=cur_starty;y<cur_endy;y++){
					pixels[cur_startx][y]->leave(x1, TOP,i);
					pixels[cur_startx][y+1]->enter(x1, BOTTOM,i);
				}
			}else { //border[bottom] down
				for(int y=cur_starty;y>cur_endy;y--){
					pixels[cur_startx][y]->leave(x1, BOTTOM,i);
					pixels[cur_startx][y-1]->enter(x1, TOP,i);
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
							pixels[x++][y]->leave(yval,RIGHT,i);
							pixels[x][y]->enter(yval,LEFT,i);
						}else{//right to left
							pixels[x--][y]->leave(yval,LEFT,i);
							pixels[x][y]->enter(yval,RIGHT,i);
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
							pixels[x][y++]->leave(xval, TOP,i);
							pixels[x][y]->enter(xval, BOTTOM,i);
						}else{// top down
							pixels[x][y--]->leave(xval, BOTTOM,i);
							pixels[x][y]->enter(xval, TOP,i);
						}
					}
				}
				// for debugging, should never happen
				if(!passed){
					vs->print();
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

	for(vector<Pixel *> &rows:pixels){
		for(Pixel *p:rows){
			for(int i=0;i<4;i++){
				if(p->intersection_nodes[i].size()>0){
					p->status = BORDER;
					break;
				}
			}
		}
	}


	for(vector<Pixel *> &rows:pixels){
		for(Pixel *pix:rows){
			pix->process_crosses(vs->num_vertices);
		}
	}
}

void MyRaster::scanline_reandering(){
	for(int y=1;y<dimy;y++){
		bool isin = false;
		for(int x=0;x<dimx;x++){
			if(pixels[x][y]->status!=BORDER){
				if(isin){
					pixels[x][y]->status = IN;
				}
				continue;
			}
			if(pixels[x][y]->intersection_nodes[BOTTOM].size()%2==1){
				isin = !isin;
			}
		}
	}
}

void MyRaster::rasterization(){

	//1. create space for the pixels
	init_pixels();

	//2. edge crossing to identify BORDER pixels
	evaluate_edges();

	//3. determine the status of rest pixels with scanline rendering
	scanline_reandering();
}

int MyRaster::get_offset_x(double xval){
	assert(mbr);
	int x = double_to_int((xval-mbr->low[0])/step_x);
	return x;
}
int MyRaster::get_offset_y(double yval){
	assert(mbr);
	int y = double_to_int((yval-mbr->low[1])/step_y);
	return y;
}
Pixel *MyRaster::get_pixel(Point &p){
	int xoff = get_offset_x(p.x);
	int yoff = get_offset_y(p.y);
	assert(xoff<=dimx);
	assert(yoff<=dimy);
	return pixels[xoff][yoff];
}
Pixel *MyRaster::get_closest_pixel(Point &p){
	int pixx = get_offset_x(p.x);
	int pixy = get_offset_y(p.y);
	if(pixx < 0){
		pixx = 0;
	}
	if(pixx > dimx){
		pixx = dimx;
	}
	if(pixy < 0){
		pixy = 0;
	}
	if(pixy > dimy){
		pixy = dimy;
	}
	return pixels[pixx][pixy];
}


vector<Pixel *> MyRaster::get_pixels(Pixel *b){

	// test all the pixels
	int txstart = get_offset_x(b->low[0]);
	int tystart = get_offset_y(b->low[1]);

	double width_d = (b->high[0]-b->low[0]+this->step_x*0.9999999)/this->step_x;
	int width = double_to_int(width_d);

	double height_d = (b->high[1]-b->low[1]+this->step_y*0.9999999)/this->step_y;
	int height = double_to_int(height_d);

	vector<Pixel *> ret;
	for(int i=txstart;i<txstart+width;i++){
		for(int j=tystart;j<tystart+height;j++){
			ret.push_back(pixels[i][j]);
		}
	}
	return ret;
}

int MyRaster::count_intersection_nodes(Point &p){
	// here we assume the point inside one of the pixel
	Pixel *pix = get_pixel(p);
	assert(pix->status==BORDER);
	int count = 0;
	for(int i=0;i<=pix->id[1];i++){
		for(double &node:pixels[pix->id[0]][i]->intersection_nodes[RIGHT]){
			if(node<=p.y){
				count++;
			}
		}
	}
	return count;
}

int MyRaster::get_num_border_edge(){
	int num = 0;
	for(vector<Pixel *> &rows:pixels){
		for(Pixel *p:rows){
			num += p->num_edges_covered();
		}
	}
	return num;
}


void MyRaster::print(){
	MyMultiPolygon *inpolys = new MyMultiPolygon();
	MyMultiPolygon *borderpolys = new MyMultiPolygon();
	MyMultiPolygon *outpolys = new MyMultiPolygon();

	for(int i=0;i<=dimx;i++){
		for(int j=0;j<=dimy;j++){
			MyPolygon *m = MyPolygon::gen_box(*pixels[i][j]);
			if(pixels[i][j]->status==BORDER){
				borderpolys->insert_polygon(m);
			}else if(pixels[i][j]->status==IN){
				inpolys->insert_polygon(m);
			}else if(pixels[i][j]->status==OUT){
				outpolys->insert_polygon(m);
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

/*
 * with the given center, expand to get the Maximum Enclosed Rectangle
 *
 * */
Pixel *MyRaster::extractMER(Pixel *starter){
	assert(starter->status==IN);
	int cx = starter->id[0];
	int cy = starter->id[1];
	Pixel *curmer = new Pixel();
	int shift[4] = {0,0,0,0};
	bool limit[4] = {false,false,false,false};
	int index = 0;
	while(!limit[0]||!limit[1]||!limit[2]||!limit[3]){
		//left
		if(!limit[0]){
			shift[0]++;
			if(cx-shift[0]<0){
				limit[0] = true;
				shift[0]--;
			}else{
				for(int i=cy-shift[1];i<=cy+shift[3];i++){
					if(pixels[cx-shift[0]][i]->status!=IN){
						limit[0] = true;
						shift[0]--;
						break;
					}
				}
			}
		}
		//bottom
		if(!limit[1]){
			shift[1]++;
			if(cy-shift[1]<0){
				limit[1] = true;
				shift[1]--;
			}else{
				for(int i=cx-shift[0];i<=cx+shift[2];i++){
					if(pixels[i][cy-shift[1]]->status!=IN){
						limit[1] = true;
						shift[1]--;
						break;
					}
				}
			}
		}
		//right
		if(!limit[2]){
			shift[2]++;
			if(cx+shift[2]>=pixels.size()){
				limit[2] = true;
				shift[2]--;
			}else{
				for(int i=cy-shift[1];i<=cy+shift[3];i++){
					if(pixels[cx+shift[2]][i]->status!=IN){
						limit[2] = true;
						shift[2]--;
						break;
					}
				}
			}
		}
		//top
		if(!limit[3]){
			shift[3]++;
			if(cy+shift[3]>=pixels[0].size()){
				limit[3] = true;
				shift[3]--;
			}else{
				for(int i=cx-shift[0];i<=cx+shift[2];i++){
					if(pixels[i][cy+shift[3]]->status!=IN){
						limit[3] = true;
						shift[3]--;
						break;
					}
				}
			}
		}
	}

	curmer->low[0] = pixels[cx-shift[0]][cy-shift[1]]->low[0];
	curmer->low[1] = pixels[cx-shift[0]][cy-shift[1]]->low[1];
	curmer->high[0] = pixels[cx+shift[2]][cy+shift[3]]->high[0];
	curmer->high[1] = pixels[cx+shift[2]][cy+shift[3]]->high[1];

	return curmer;
}


bool MyRaster::contain(Pixel *b, bool &contained){

	if(!mbr->contain(*b)){
		contained = false;
		return true;
	}
	// test all the pixels
	vector<Pixel *> covered = get_pixels(b);

	int incount = 0;
	int outcount = 0;
	for(Pixel *pix:covered){
		if(pix->status==OUT){
			outcount++;
		}
		if(pix->status==IN){
			incount++;
		}
	}
	int total = covered.size();
	covered.clear();
	// all in/out
	if(incount==total){
		contained = true;
		return true;
	}else if(outcount==total){
		contained = false;
		return true;
	}

	// not determined by checking pixels only
	return false;
}



int MyRaster::get_num_pixels(PartitionStatus status){
	int num = 0;
	for(vector<Pixel *> &rows:pixels){
		for(Pixel *p:rows){
			if(p->status==status){
				num++;
			}
		}
	}
	return num;
}

size_t MyRaster::get_num_pixels(){
	return dimx*dimy;
}

size_t MyRaster::get_num_gridlines(){
	return dimx+dimy;
}

size_t MyRaster::get_num_crosses(){
	size_t num = 0;
	for(vector<Pixel *> &rows:pixels){
		for(Pixel *p:rows){
			num += p->intersection_nodes[RIGHT].size();
			num += p->intersection_nodes[BOTTOM].size();
			num += p->intersection_nodes[LEFT].size();
			num += p->intersection_nodes[TOP].size();
		}
	}
	return num;
}

vector<Pixel *> MyRaster::get_pixels(PartitionStatus status){
	vector<Pixel *> ret;
	for(vector<Pixel *> &rows:pixels){
		for(Pixel *p:rows){
			if(p->status==status){
				ret.push_back(p);
			}
		}
	}
	return ret;
}
