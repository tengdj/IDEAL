/*
 * Pixel.cpp
 *
 *  Created on: May 26, 2020
 *      Author: teng
 */


#include "Pixel.h"

bool print_debug = false;

void box::update(Point &p){
	if(low[0]>p.x){
		low[0] = p.x;
	}
	if(high[0]<p.x){
		high[0] = p.x;
	}

	if(low[1]>p.y){
		low[1] = p.y;
	}
	if(high[1]<p.y){
		high[1] = p.y;
	}
}

bool box::intersect(box &target){
	return !(target.low[0]>high[0]||
			 target.high[0]<low[0]||
			 target.low[1]>high[1]||
			 target.high[1]<low[1]);
}
bool box::contain(box &target){
	return target.low[0]>=low[0]&&
		   target.high[0]<=high[0]&&
		   target.low[1]>=low[1]&&
		   target.high[1]<=high[1];
}
bool box::contain(Point &p){
	return p.x>=low[0]&&
		   p.x<=high[0]&&
		   p.y>=low[1]&&
		   p.y<=high[1];
}

double box::max_distance(Point &p){
	double md = 0;
	double dist = (p.x-low[0])*(p.x-low[0])+(p.y-low[1])*(p.y-low[1]);
	if(dist>md){
		md = dist;
	}
	dist = (p.x-low[0])*(p.x-low[0])+(p.y-high[1])*(p.y-high[1]);
	if(dist>md){
		md = dist;
	}
	dist = (p.x-high[0])*(p.x-high[0])+(p.y-low[1])*(p.y-low[1]);
	if(dist>md){
		md = dist;
	}
	dist = (p.x-high[0])*(p.x-high[0])+(p.y-high[1])*(p.y-high[1]);
	if(dist>md){
		md = dist;
	}
	return sqrt(md);
}

double box::area(){
	return (high[0]-low[0])*(high[1]-low[1]);
}



double box::distance(Point &p){
	if(this->contain(p)){
		return 0;
	}
	double dx = max(abs(p.x-(low[0]+high[0])/2) - (high[0]-low[0])/2, 0.0);
	double dy = max(abs(p.y-(low[1]+high[1])/2) - (high[1]-low[1])/2, 0.0);
	return sqrt(dx * dx + dy * dy);
}

double box::distance_geography(Point &p){
	if(this->contain(p)){
		return 0.0;
	}
	double dx = max(abs(p.x-(low[0]+high[0])/2) - (high[0]-low[0])/2, 0.0);
	double dy = max(abs(p.y-(low[1]+high[1])/2) - (high[1]-low[1])/2, 0.0);
	dx = dx/degree_per_kilometer_latitude;
	dy = dy/degree_per_kilometer_longitude(p.y);
	return sqrt(dx * dx + dy * dy);
}

void box::print_vertices(){
	printf("%f %f, %f %f, %f %f, %f %f, %f %f",
				low[0],low[1],
				high[0],low[1],
				high[0],high[1],
				low[0],high[1],
				low[0],low[1]);
}

void box::print(){
	printf("POLYGON((");
	print_vertices();
	printf("))\n");

}


void Pixel::enter(double val, Direction d, int vnum){
	cross_info ci;
	ci.type = ENTER;
	ci.vertex = val;
	crosses[d].push_back(ci);
	if(vstart==-1){
		vstart = vnum;
	}
	vend = vnum;
}



void Pixel::leave(double val, Direction d, int vnum){
	cross_info ci;
	ci.type = LEAVE;
	ci.vertex = val;
	crosses[d].push_back(ci);
	if(vstart==-1){
		vstart = vnum;
	}
	vend = vnum;
}

