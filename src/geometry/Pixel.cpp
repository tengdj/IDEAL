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

/*
 * distance related functions
 * */
double box::distance(box &t, bool geography){
	if(this->intersect(t)){
		return 0;
	}
	double dx = 0;
	double dy = 0;

	if(t.low[1]>high[1]){
		dy = t.low[1] - high[1];
	}else if(t.high[1]<low[1]){
		dy = low[1] - t.high[1];
	}else{
		dy = 0;
	}

	if(t.low[0]>high[0]){
		dx = t.low[0] - high[0];
	}else if(t.high[0]<low[0]){
		dx = low[0] - t.high[0];
	}else{
		dx = 0;
	}

	if(geography){
		dy = dy/degree_per_kilometer_latitude;
		dx = dx/degree_per_kilometer_longitude(low[1]);
	}
	return sqrt(dx * dx + dy * dy);
}

double box::max_distance(box &target, bool geography){
	double highx = max(high[0],target.high[0]);
	double highy = max(high[1],target.high[1]);
	double lowx = min(low[0], target.low[0]);
	double lowy = min(low[1], target.low[1]);
	double dx = highx - lowx;
	double dy = highy - lowy;
	if(geography){
		dy = dy/degree_per_kilometer_latitude;
		dx = dx/degree_per_kilometer_longitude(low[1]);
	}
	return sqrt(dx * dx + dy * dy);
}

double box::distance(Point &p, bool geography){
	if(this->contain(p)){
		return 0;
	}
	double dx = max(abs(p.x-(low[0]+high[0])/2) - (high[0]-low[0])/2, 0.0);
	double dy = max(abs(p.y-(low[1]+high[1])/2) - (high[1]-low[1])/2, 0.0);
	if(geography){
		dy = dy/degree_per_kilometer_latitude;
		dx = dx/degree_per_kilometer_longitude(p.y);
	}
	return sqrt(dx * dx + dy * dy);
}

double box::max_distance(Point &p, bool geography){

	double dx = max(abs(p.x-low[0]), abs(p.x-high[0]));
	double dy = max(abs(p.y-low[1]), abs(p.y-high[1]));

	if(geography){
		dy = dy/degree_per_kilometer_latitude;
		dx = dx/degree_per_kilometer_longitude(p.y);
	}
	return sqrt(dx*dx+dy*dy);
}

/*
 * print functions
 * */
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

	//very very very very rare cases
	if(crosses.size()%2==1){
		crosses.push_back(cross_info((cross_type)!crosses[crosses.size()-1].type,crosses[crosses.size()-1].edge_id));
	}

	assert(crosses.size()%2==0);
	int start = 0;
	int end = crosses.size()-1;

	//special case for the first edge
	if(crosses[0].type==LEAVE){
		assert(crosses[end].type==ENTER);
		edge_ranges.push_back(edge_range(crosses[end].edge_id,num_edges-2));
		edge_ranges.push_back(edge_range(0,crosses[0].edge_id));
		start++;
		end--;
	}

	for(int i=start;i<=end;i++){
		assert(crosses[i].type==ENTER);
		//special case, an ENTER has no pair LEAVE,
		//happens when one edge crosses the pair
		if(i==end||crosses[i+1].type==ENTER){
			edge_ranges.push_back(edge_range(crosses[i].edge_id,crosses[i].edge_id));
		}else{
			edge_ranges.push_back(edge_range(crosses[i].edge_id,crosses[i+1].edge_id));
			i++;
		}
	}

	// confirm the correctness
	for(edge_range &r:edge_ranges){
		assert(r.vstart<=r.vend&&r.vend<num_edges);
	}
	crosses.clear();
}

int Pixel::num_edges_covered(){
	int c = 0;
	for(edge_range &r:edge_ranges){
		c += r.vend-r.vstart+1;
	}
	return c;
}
