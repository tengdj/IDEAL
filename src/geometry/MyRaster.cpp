#include "../include/MyRaster.h"
#include "../include/MyPolygon.h"

void MyRaster::init_raster(int num_pixels){
    assert(num_pixels > 0);

	double multi = abs((mbr->high[1]-mbr->low[1])/(mbr->high[0]-mbr->low[0]));
	dimx = std::pow((num_pixels)/multi,0.5);
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