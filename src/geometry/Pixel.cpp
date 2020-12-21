/*
 * Pixel.cpp
 *
 *  Created on: May 26, 2020
 *      Author: teng
 */


#include "MyPolygon.h"

bool print_debug = false;




void Pixel::enter(double val, Direction d, int vnum){
	if(print_debug){
		cout<<direction_str[d];
		cout<<" enter "<<id[0]<<" "<<id[1]<<endl;
	}
	cross_info ci;
	ci.type = ENTER;
	ci.vertex = val;
	ci.direction = d;
	crosses.push_back(ci);
	if(vstart==-1){
		vstart = vnum;
	}
	vend = vnum;
}



void Pixel::leave(double val, Direction d, int vnum){
	if(print_debug){
		cout<<direction_str[d];
		cout<<" leave "<<id[0]<<" "<<id[1]<<endl;
	}
	cross_info ci;
	ci.type = LEAVE;
	ci.vertex = val;
	ci.direction = d;
	crosses.push_back(ci);
	if(vstart==-1){
		vstart = vnum;
	}
	vend = vnum;
}


MyPolygon *Pixel::to_polygon(){
	return 	MyPolygon::gen_box(low[0],low[1],high[0],high[1]);
};

