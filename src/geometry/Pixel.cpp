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

inline double geogdist(double x1, double y1, double x2, double y2){

	double dx = x2-x1;
	double dy = y2-y1;
	dx = dx/degree_per_kilometer_latitude;
	dy = dy/degree_per_kilometer_longitude(y1);
	return dx*dx+dy*dy;
}

double box::max_distance_geography(Point &p){
	double md = 0;
	double dist = geogdist(p.x, p.y, low[0], low[1]);
	md = dist;
	dist = geogdist(p.x, p.y, low[0], high[1]);
	if(dist>md){
		md = dist;
	}
	dist = geogdist(p.x, p.y, high[0], low[1]);
	if(dist>md){
		md = dist;
	}
	dist = geogdist(p.x, p.y, high[0], high[1]);
	if(dist>md){
		md = dist;
	}
	return sqrt(md);
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
	intersection_nodes[d].push_back(val);
	crosses.push_back(cross_info(ENTER,vnum));
}



void Pixel::leave(double val, Direction d, int vnum){
	intersection_nodes[d].push_back(val);
	crosses.push_back(cross_info(LEAVE,vnum));
}

void Pixel::process_crosses(int num_edges){
	if(crosses.size()==0){
		return;
	}
	vstart = 0;
	vend = num_edges-1;

	vstart = crosses[0].edge_num;
	vend = crosses[crosses.size()-1].edge_num;

	//edge_ranges.push_back(edge_range(vstart,vend));return;
	//very very very very rare cases
	if(crosses.size()%2==1){
		crosses.push_back(cross_info((cross_type)!crosses[crosses.size()-1].type,crosses[crosses.size()-1].edge_num));
	}

	assert(crosses.size()%2==0);
	int start = 0;
	int end = crosses.size()-1;

	//special case
	if(crosses[0].type==LEAVE){
		assert(crosses[end].type==ENTER);
		edge_ranges.push_back(edge_range(crosses[end].edge_num,num_edges-2));
		edge_ranges.push_back(edge_range(0,crosses[0].edge_num));
		start++;
		end--;
	}

	for(int i=start;i<=end;i++){
		assert(crosses[i].type==ENTER);
		//special case, an ENTER has no pair LEAVE
		if(i==end||crosses[i+1].type==ENTER){
			edge_ranges.push_back(edge_range(crosses[i].edge_num,crosses[i].edge_num));
		}else{
			edge_ranges.push_back(edge_range(crosses[i].edge_num,crosses[i+1].edge_num));
			i++;
		}
	}

	for(edge_range &r:edge_ranges){
		assert(r.vstart<=r.vend&&r.vend<num_edges);
	}
}

int Pixel::num_edges_covered(){
		int c = 0;
		for(edge_range &r:edge_ranges){
			c += r.vend-r.vstart+1;
		}
		return c;
	}


