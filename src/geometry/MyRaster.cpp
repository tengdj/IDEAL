/*
 * MyRaster.cpp
 *
 *  Created on: Jan 4, 2021
 *      Author: teng
 */


#include "../include/MyPolygon.h"

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
	if(mbr != nullptr) delete mbr;
	if(pixs != nullptr) delete pixs;
	if(horizontal != nullptr) delete horizontal;
	if(vertical != nullptr) delete vertical;
}

void MyRaster::process_crosses(map<int, vector<cross_info>> edges_info){
	int num_edge_seqs = 0;
	for(auto ei : edges_info){
		num_edge_seqs += ei.second.size();
	}
	pixs->init_edge_sequences(num_edge_seqs);

	int idx = 0;
	int edge_count = 0;
	for(auto info : edges_info){
		auto pix = info.first;
		auto crosses = info.second; 
		if(crosses.size() == 0) return;
		
		if(crosses.size() % 2 == 1){
			crosses.push_back(cross_info((cross_type)!crosses[crosses.size()-1].type, crosses[crosses.size()-1].edge_id));
		}
		
		assert(crosses.size()%2==0);

		// Initialize based on crosses.size().
		int start = 0;
		int end = crosses.size() - 1;
		pixs->set_offset(pix, idx);

		if(crosses[0].type == LEAVE){
			assert(crosses[end].type == ENTER);
			pixs->add_edge(idx ++, 0, crosses[0].edge_id);
			pixs->add_edge(idx ++, crosses[end].edge_id, vs->num_vertices - 2);
			start ++;
			end --;
		}

		for(int i = start; i <= end; i++){
			assert(crosses[i].type == ENTER);
			//special case, an ENTER has no pair LEAVE,
			//happens when one edge crosses the pair
			if(i == end || crosses[i + 1].type == ENTER){
				pixs->add_edge(idx ++, crosses[i].edge_id, crosses[i].edge_id);
			}else{
				pixs->add_edge(idx ++, crosses[i].edge_id, crosses[i+1].edge_id);
				i++;
			}
		}
	}
}

void MyRaster::process_intersection(map<int, vector<double>> intersection_info, string direction){
	int num_nodes = 0;
	for(auto i : intersection_info){
		num_nodes += i.second.size();
	}
	if(direction == "horizontal"){
		horizontal->init_intersection_node(num_nodes);
		horizontal->set_num_crosses(num_nodes);
		int idx = 0;
		for(auto info : intersection_info){
			auto h = info.first;
			auto nodes = info.second;
			
			sort(nodes.begin(), nodes.end());

			horizontal->set_offset(h, idx);

			for(auto node : nodes){
				horizontal->add_node(idx, node);
				idx ++;
			}
		}
		horizontal->set_offset(dimy, idx);
	}else{
		vertical->init_intersection_node(num_nodes);
		vertical->set_num_crosses(num_nodes);

		int idx = 0;
		for(auto info : intersection_info){
			auto h = info.first;
			auto nodes = info.second;
			
			sort(nodes.begin(), nodes.end());

			vertical->set_offset(h, idx);

			for(auto node : nodes){
				vertical->add_node(idx, node);
				idx ++;
			}
		}
		vertical->set_offset(dimx, idx);		
	}

}

void MyRaster::init_pixels(){
	assert(mbr);
	pixs = new Pixels((dimx+1)*(dimy+1));
	horizontal = new Grid_line(dimy);
	vertical = new Grid_line(dimx);
}

void MyRaster::evaluate_edges(){
	map<int, vector<double>> horizontal_intersect_info;
	map<int, vector<double>> vertical_intersect_info;
	map<int, vector<cross_info>> edges_info;
	
	// normalize
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
		// todo should not happen for normal cases
		if(cur_startx>dimx||cur_endx>dimx||cur_starty>dimy||cur_endy>dimy){
			cout<<"xrange\t"<<cur_startx<<" "<<cur_endx<<endl;
			cout<<"yrange\t"<<cur_starty<<" "<<cur_endy<<endl;
			printf("xrange_val\t%f %f\n",(x1-start_x)/step_x, (x2-start_x)/step_x);
			printf("yrange_val\t%f %f\n",(y1-start_y)/step_y, (y2-start_y)/step_y);
			assert(false);
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
					vertical_intersect_info[x + 1].push_back(y1);
					edges_info[get_id(x, cur_starty)].push_back(cross_info(LEAVE, i));
					edges_info[get_id(x + 1, cur_starty)].push_back(cross_info(ENTER, i));
					pixs->set_status(get_id(x, cur_starty), BORDER);
					pixs->set_status(get_id(x+1, cur_starty), BORDER);
				}
			}else { // right to left
				for(int x=cur_startx;x>cur_endx;x--){
					vertical_intersect_info[x].push_back(y1);
					edges_info[get_id(x, cur_starty)].push_back(cross_info(LEAVE, i));
					edges_info[get_id(x - 1, cur_starty)].push_back(cross_info(ENTER, i));
					pixs->set_status(get_id(x, cur_starty), BORDER);
					pixs->set_status(get_id(x-1, cur_starty), BORDER);
				}
			}
		}else if(x1==x2){
			//bottom up
			if(cur_starty<cur_endy){
				for(int y=cur_starty;y<cur_endy;y++){
					horizontal_intersect_info[y + 1].push_back(x1);
					edges_info[get_id(cur_startx, y)].push_back(cross_info(LEAVE, i));
					edges_info[get_id(cur_startx, y + 1)].push_back(cross_info(ENTER, i));
					pixs->set_status(get_id(cur_startx, y), BORDER);
					pixs->set_status(get_id(cur_startx, y+1), BORDER);
				}
			}else { //border[bottom] down
				for(int y=cur_starty;y>cur_endy;y--){
					horizontal_intersect_info[y].push_back(x1);
					edges_info[get_id(cur_startx, y)].push_back(cross_info(LEAVE, i));
					edges_info[get_id(cur_startx, y - 1)].push_back(cross_info(ENTER, i));
					pixs->set_status(get_id(cur_startx, y), BORDER);
					pixs->set_status(get_id(cur_startx, y-1), BORDER);
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
							vertical_intersect_info[x + 1].push_back(yval);
							pixs->set_status(get_id(x, y), BORDER);
							edges_info[get_id(x ++, y)].push_back(cross_info(LEAVE, i));
							edges_info[get_id(x, y)].push_back(cross_info(ENTER, i));
							pixs->set_status(get_id(x, y), BORDER);
						}else{//right to left
							vertical_intersect_info[x].push_back(yval);
							pixs->set_status(get_id(x, y), BORDER);
							edges_info[get_id(x --, y)].push_back(cross_info(LEAVE, i));
							edges_info[get_id(x, y)].push_back(cross_info(ENTER, i));
							pixs->set_status(get_id(x, y), BORDER);
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
							horizontal_intersect_info[y + 1].push_back(xval);
							pixs->set_status(get_id(x, y), BORDER);
							edges_info[get_id(x, y ++)].push_back(cross_info(LEAVE, i));
							edges_info[get_id(x, y)].push_back(cross_info(ENTER, i));
							pixs->set_status(get_id(x, y), BORDER);
						}else{// top down
							horizontal_intersect_info[y].push_back(xval);
							pixs->set_status(get_id(x, y), BORDER);
							edges_info[get_id(x, y --)].push_back(cross_info(LEAVE, i));
							edges_info[get_id(x, y)].push_back(cross_info(ENTER, i));
							pixs->set_status(get_id(x, y), BORDER);
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


	process_crosses(edges_info);
	process_intersection(horizontal_intersect_info, "horizontal");
	process_intersection(vertical_intersect_info, "vertical");
	pixs->process_pixels_null(dimx, dimy);
}

void MyRaster::scanline_reandering(){
	const double start_x = mbr->low[0];
	const double start_y = mbr->low[1];

	for(int y = 1; y < dimy; y ++){
		bool isin = false;
		uint16_t i = horizontal->get_offset(y), j = horizontal->get_offset(y + 1);
		for(int x = 0; x < dimx; x ++){
			if(pixs->show_status(get_id(x, y)) != BORDER){
				if(isin){
					pixs->set_status(get_id(x, y), IN);
				}else{
					pixs->set_status(get_id(x, y), OUT);
				}
				continue;
			}
			int pass = 0;
			while(i < j && horizontal->get_intersection_nodes(i) <= start_x + step_x * (x + 1)){
				pass ++;
				i ++;
			}
			if(pass % 2 == 1) isin = !isin;

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

uint16_t MyRaster::get_num_sequences(int id){
	if(pixs->show_status(id) != BORDER) return 0;
	return pixs->get_pointer(id + 1) - pixs->get_pointer(id);
}

// similar to the expand_radius function, get the possible minimum distance between point p and
// the target pixels which will be evaluated in this step, will be used as a stop sign
double MyRaster::get_possible_min(Point &p, int center, int step, bool geography){
	int core_x_low = get_x(center);
	int core_x_high = get_x(center);
	int core_y_low = get_y(center);
	int core_y_high = get_y(center);

	vector<int> needprocess;

	int ymin = max(0,core_y_low-step);
	int ymax = min(dimy,core_y_high+step);

	double mindist = DBL_MAX;
	//left scan
	if(core_x_low-step>=0){
		double x = get_pixel_box(core_x_low-step,ymin).high[0];
		double y1 = get_pixel_box(core_x_low-step,ymin).low[1];
		double y2 = get_pixel_box(core_x_low-step,ymax).high[1];

		Point p1 = Point(x, y1);
		Point p2 = Point(x, y2);
		double dist = point_to_segment_distance(p, p1, p2, geography);
		mindist = min(dist, mindist);
	}
	//right scan
	if(core_x_high+step<=get_dimx()){
		double x = get_pixel_box(core_x_high+step,ymin).low[0];
		double y1 = get_pixel_box(core_x_high+step,ymin).low[1];
		double y2 = get_pixel_box(core_x_high+step,ymax).high[1];
		Point p1 = Point(x, y1);
		Point p2 = Point(x, y2);
		double dist = point_to_segment_distance(p, p1, p2, geography);
		mindist = min(dist, mindist);
	}

	// skip the first if there is left scan
	int xmin = max(0,core_x_low-step+(core_x_low-step>=0));
	// skip the last if there is right scan
	int xmax = min(dimx,core_x_high+step-(core_x_high+step<=dimx));
	//bottom scan
	if(core_y_low-step>=0){
		double y = get_pixel_box(xmin,core_y_low-step).high[1];
		double x1 = get_pixel_box(xmin,core_y_low-step).low[0];
		double x2 = get_pixel_box(xmax,core_y_low-step).high[0];
		Point p1 = Point(x1, y);
		Point p2 = Point(x2, y);
		double dist = point_to_segment_distance(p, p1, p2, geography);
		mindist = min(dist, mindist);
	}
	//top scan
	if(core_y_high+step<=get_dimy()){
		double y = get_pixel_box(xmin,core_y_low+step).low[1];
		double x1 = get_pixel_box(xmin,core_y_low+step).low[0];
		double x2 = get_pixel_box(xmax,core_y_low+step).high[0];
		Point p1 = Point(x1, y);
		Point p2 = Point(x2, y);
		double dist = point_to_segment_distance(p, p1, p2, geography);
		mindist = min(dist, mindist);
	}
	return mindist;
}

vector<int> MyRaster::expand_radius(int center, int step){

	int lowx = get_x(center);
	int highx = get_x(center);
	int lowy = get_y(center);
	int highy = get_y(center);

	return expand_radius(lowx,highx,lowy,highy,step);
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

// retrieve the pixels in the raster which is closest to the target pixels
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

// retrieve the pixel in the raster which is closest to the target pixels
// pick the one in the middle if there is multiple pixels

int MyRaster::count_intersection_nodes(Point &p){
	// here we assume the point inside one of the pixel
	int pix_id = get_pixel_id(p);
	assert(pixs->show_status(pix_id) == BORDER);
	int count = 0;
	int x = get_x(pix_id) + 1;
	uint16_t i = vertical->get_offset(x), j;
	if(x < dimx) j = vertical->get_offset(x + 1);
	else j = vertical->get_num_crosses();
	while(i < j && vertical->get_intersection_nodes(i) <= p.y){
		count ++;
		i ++;
	}
	return count;
}



void MyRaster::print(){
	MyMultiPolygon *inpolys = new MyMultiPolygon();
	MyMultiPolygon *borderpolys = new MyMultiPolygon();
	MyMultiPolygon *outpolys = new MyMultiPolygon();

	for(int i=0;i<=dimx;i++){
		for(int j=0;j<=dimy;j++){
			box bx = get_pixel_box(i, j);
			MyPolygon *m = MyPolygon::gen_box(bx);
			if(pixs->show_status(get_id(i, j)) == BORDER){
				borderpolys->insert_polygon(m);
			}else if(pixs->show_status(get_id(i, j)) == IN){
				inpolys->insert_polygon(m);
			}else if(pixs->show_status(get_id(i, j)) == OUT){
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

/*
 * with the given center, expand to get the Maximum Enclosed Rectangle
 *
 * */
box *MyRaster::extractMER(int starter){
	assert(pixs->show_status(starter) == IN);
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
					if(pixs->show_status(get_id(0, i) != IN))
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
					if(pixs->show_status(get_id(i, cy-shift[1]) != IN)){
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

size_t MyRaster::get_num_gridlines(){
	return dimx+dimy;
}