#include "../include/MyRaster.h"

void MyRaster::init_raster(int num_pixels){
    assert(num_pixels >= 0);

	double multi = abs((mbr->high[1]-mbr->low[1])/(mbr->high[0]-mbr->low[0]));
	dimx = std::pow(num_pixels/multi,0.5);
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

	status = new uint8_t[(dimx+1)*(dimy+1) / 4 + 1];
    memset(status, 0, ((dimx+1)*(dimy+1) / 4 + 1) * sizeof(uint8_t));
}

void MyRaster::rasterization(VertexSequence *vs, int vpr){
	pthread_mutex_lock(&raster_lock);
	mbr = vs->getMBR();
	init_raster(vs->num_vertices / vpr);
	rasterization(vs);
	pthread_mutex_unlock(&raster_lock);
}

void MyRaster::rasterization(VertexSequence *vs){
	assert(mbr);
	const double start_x = mbr->low[0];
	const double start_y = mbr->low[1];

	for(int i=0;i<vs->num_vertices-1;i++){
		double x1 = vs->p[i].x;
		double y1 = vs->p[i].y;
		double x2 = vs->p[i+1].x;
		double y2 = vs->p[i+1].y;

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

		assert(cur_startx<=dimx);
		assert(cur_endx<=dimx);
		assert(cur_starty<=dimy);
		assert(cur_endy<=dimy);

		//in the same pixel
		if(cur_startx==cur_endx&&cur_starty==cur_endy){
			continue;
		}

		if(y1==y2){
			//left to right
			if(cur_startx<cur_endx){
				for(int x=cur_startx;x<cur_endx;x++){
					set_status(get_id(x, cur_starty), BORDER);
					set_status(get_id(x+1, cur_starty), BORDER);
				}
			}else { // right to left
				for(int x=cur_startx;x>cur_endx;x--){
					set_status(get_id(x, cur_starty), BORDER);
					set_status(get_id(x-1, cur_starty), BORDER);
				}
			}
		}else if(x1==x2){
			//bottom up
			if(cur_starty<cur_endy){
				for(int y=cur_starty;y<cur_endy;y++){
					set_status(get_id(cur_startx, y), BORDER);
					set_status(get_id(cur_startx, y+1), BORDER);
				}
			}else { //border[bottom] down
				for(int y=cur_starty;y>cur_endy;y--){
					set_status(get_id(cur_startx, y), BORDER);
					set_status(get_id(cur_startx, y-1), BORDER);
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
							set_status(get_id(x ++, y), BORDER);
							set_status(get_id(x, y), BORDER);
						}else{//right to left
							set_status(get_id(x --, y), BORDER);
							set_status(get_id(x, y), BORDER);
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
							set_status(get_id(x, y++), BORDER);
							set_status(get_id(x, y), BORDER);
						}else{// top down
							set_status(get_id(x, y --), BORDER);
							set_status(get_id(x, y), BORDER);
						}
					}
				}
			}
		}
	}



	for(int i = 0; i < get_num_pixels(); i ++){
		if(show_status(i) == BORDER) continue;
		box bx = get_pixel_box(get_x(i), get_y(i));
		Point point;
		point.set(bx.low[0], bx.low[1]);
		if(vs->contain(point)) set_status(i, IN);
		else set_status(i, OUT);
	}
}

void MyRaster::print(){
	MyMultiPolygon *inpolys = new MyMultiPolygon();
	MyMultiPolygon *borderpolys = new MyMultiPolygon();
	MyMultiPolygon *outpolys = new MyMultiPolygon();

	for(int i=0;i<=dimx;i++){
		for(int j=0;j<=dimy;j++){
			box bx = get_pixel_box(i, j);
			MyPolygon *m = MyPolygon::gen_box(bx);
			if(show_status(get_id(i, j)) == BORDER){
				borderpolys->insert_polygon(m);
			}else if(show_status(get_id(i, j)) == IN){
				inpolys->insert_polygon(m);
			}else if(show_status(get_id(i, j)) == OUT){
				outpolys->insert_polygon(m);
			}
		}
	}

	cout<<"border:" << borderpolys->num_polygons() <<endl;
	borderpolys->print();
	cout<<"in:"<< inpolys->num_polygons() << endl;
	inpolys->print();
	cout<<"out:"<< outpolys->num_polygons() << endl;
	outpolys->print();
	cout << endl;


	delete borderpolys;
	delete inpolys;
	delete outpolys;
}

box *MyRaster::extractMER(int starter){
	assert(show_status(starter) == IN);
	int cx = get_x(starter);
	int cy = get_y(starter);
	box *curmer = new box();
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
					//log("%d %d %d %d", cx-shift[0], dimx, i, dimy);
					if(show_status(get_id(cx-shift[0], i)) != IN){
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
					if(show_status(get_id(i, cy-shift[1])) != IN){
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
			if(cx+shift[2]>=(dimx+1)){
				limit[2] = true;
				shift[2]--;
			}else{
				for(int i=cy-shift[1];i<=cy+shift[3];i++){
					if(show_status(get_id(cx+shift[2], i)) != IN){
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
			if(cy+shift[3]>=(dimy+1)){
				limit[3] = true;
				shift[3]--;
			}else{
				for(int i=cx-shift[0];i<=cx+shift[2];i++){
					// log("%d %d %d %d", i, dimx, cy+shift[3], dimy);
					if(show_status(get_id(i, cy+shift[3])) != IN){
						limit[3] = true;
						shift[3]--;
						break;
					}
				}
			}
		}
	}

	curmer->low[0] = get_pixel_box(cx-shift[0], cy-shift[1]).low[0];
	curmer->low[1] = get_pixel_box(cx-shift[0], cy-shift[1]).low[1];
	curmer->high[0] = get_pixel_box(cx+shift[2], cy+shift[3]).high[0];
	curmer->high[1] = get_pixel_box(cx+shift[2], cy+shift[3]).high[1];

	return curmer;
}

int MyRaster::get_id(int x, int y){
	assert(x>=0&&x<=dimx);
	assert(y>=0&&y<=dimy);
	return y * (dimx+1) + x;
}

// from id to pixel x
int MyRaster::get_x(int id){
	return id % (dimx+1);
}

// from id to pixel y
int MyRaster::get_y(int id){
	assert((id / (dimx+1)) <= dimy);
	return id / (dimx+1);
}

// the range must be [0, dimx]
int MyRaster::get_offset_x(double xval){
	assert(mbr);
	assert(step_x>0.000000000001 && "the width per pixel must not be 0");
	int x = double_to_int((xval-mbr->low[0])/step_x);
	return min(max(x, 0), dimx);
}
// the range must be [0, dimy]
int MyRaster::get_offset_y(double yval){
	assert(mbr);
	assert(step_y>0.000000000001 && "the hight per pixel must not be 0");
	int y = double_to_int((yval-mbr->low[1])/step_y);
	return min(max(y, 0), dimy);
}

void MyRaster::set_status(int id, PartitionStatus state){
	int pos = id % 4 * 2;   // The multiplication by 2 is because each status occupies 2 bits.
	if(state == OUT){
		status[id / 4] &= ~((uint8_t)3 << pos);
	}else if(state == IN){
		status[id / 4] |= ((uint8_t)3 << pos);
	}else{
		status[id / 4] &= ~((uint8_t)1 << pos);
		status[id / 4] |= ((uint8_t)1 << (pos + 1));
	}
}

PartitionStatus MyRaster::show_status(int id){
	uint8_t st = status[id / 4];
	int pos = id % 4 * 2;   // The multiplication by 2 is because each status occupies 2 bits.	
	st &= ((uint8_t)3 << pos);
	st >>= pos;
	if(st == 0) return OUT;
	if(st == 3) return IN;
	return BORDER;
}

vector<int> MyRaster::get_closest_pixels(box &target){

	// note that at 0 or dimx/dimy will be returned if
	// the range of target is beyound this, as expected
	int txstart = get_offset_x(target.low[0]);
	int txend = get_offset_x(target.high[0]);
	int tystart = get_offset_y(target.low[1]);
	int tyend = get_offset_y(target.high[1]);
	vector<int> ret;
	for(int i=txstart;i<=txend;i++){
		for(int j=tystart;j<=tyend;j++){
			ret.push_back(get_id(i, j));
		}
	}
	return ret;
}

int MyRaster::get_closest_pixel(Point &p){
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
	return get_id(pixx, pixy);
}

vector<int> MyRaster::get_pixels(PartitionStatus status){
	vector<int> ret;
	for(int id = 0; id < get_num_pixels(); id ++){
		if(show_status(id) == status){
			ret.push_back(id);
		}
	}
	return ret;
}

box MyRaster::get_pixel_box(int x, int y){
	const double start_x = mbr->low[0];
	const double start_y = mbr->low[1];

	double lowx = start_x + x * step_x;
	double lowy = start_y + y * step_y;
	double highx = start_x + (x + 1) * step_x;
	double highy = start_y + (y + 1) * step_y;

	return box(lowx, lowy, highx, highy);
}

int MyRaster::get_pixel_id(Point &p){
	int xoff = get_offset_x(p.x);
	int yoff = get_offset_y(p.y);
	assert(xoff <= dimx);
	assert(yoff <= dimy);
	return get_id(xoff, yoff);
}

vector<int> MyRaster::retrieve_pixels(box *target){
	vector<int> ret;
	int start_x = get_offset_x(target->low[0]);
	int start_y = get_offset_y(target->low[1]);
	int end_x = get_offset_x(target->high[0]);
	int end_y = get_offset_y(target->high[1]);

	//log("%d %d %d %d %d %d",dimx,dimy,start_x,end_x,start_y,end_y);
	for(int i=start_x;i<=end_x;i++){
		for(int j=start_y;j<=end_y;j++){
			ret.push_back(get_id(i , j));
		}
	}
	return ret;
}

vector<int> MyRaster::expand_radius(int core_x_low, int core_x_high, int core_y_low, int core_y_high, int step){

	vector<int> needprocess;
	int ymin = max(0,core_y_low-step);
	int ymax = min(dimy,core_y_high+step);

	//left scan
	if(core_x_low-step>=0){
		for(int y=ymin;y<=ymax;y++){
			needprocess.push_back(get_id(core_x_low-step,y));
		}
	}
	//right scan
	if(core_x_high+step<=get_dimx()){
		for(int y=ymin;y<=ymax;y++){
			needprocess.push_back(get_id(core_x_high+step,y));
		}
	}

	// skip the first if there is left scan
	int xmin = max(0,core_x_low-step+(core_x_low-step>=0));
	// skip the last if there is right scan
	int xmax = min(dimx,core_x_high+step-(core_x_high+step<=dimx));
	//bottom scan
	if(core_y_low-step>=0){
		for(int x=xmin;x<=xmax;x++){
			needprocess.push_back(get_id(x,core_y_low-step));
		}
	}
	//top scan
	if(core_y_high+step<=get_dimy()){
		for(int x=xmin;x<=xmax;x++){
			needprocess.push_back(get_id(x,core_y_high+step));
		}
	}

	return needprocess;
}

vector<int> MyRaster::expand_radius(int center, int step){

	int lowx = get_x(center);
	int highx = get_x(center);
	int lowy = get_y(center);
	int highy = get_y(center);

	return expand_radius(lowx,highx,lowy,highy,step);
}

size_t MyRaster::get_num_pixels(){
	return (dimx+1)*(dimy+1);
}

size_t MyRaster::get_num_pixels(PartitionStatus status){
	size_t num = 0;
	for(int i = 0; i < get_num_pixels();i ++){
		if(show_status(i) == status){
			num++;
		}
	}
	return num;	
}