/*
 * vertex_sequence.cpp
 *
 *  Created on: Apr 25, 2023
 *      Author: teng
 */




#include "../include/MyPolygon.h"

VertexSequence::VertexSequence(int nv){
	p = new Point[nv];
	num_vertices = nv;
}
VertexSequence::VertexSequence(int nv, Point *pp){
	assert(nv>0);
	num_vertices = nv;
	p = new Point[nv];
	memcpy((char *)p,(char *)pp,num_vertices*sizeof(Point));
};

VertexSequence::~VertexSequence(){
	if(p){
		delete []p;
	}
}

vector<Vertex *> VertexSequence::pack_to_polyline(){
	vector<Vertex *> polyline;
	polyline.resize(num_vertices-1);
	for(int i=0;i<num_vertices-1;i++){
		polyline[i] = new Vertex(p[i].x,p[i].y);
	}
	return polyline;
}

box *VertexSequence::getMBR(){
	box *mbr = new box();
	for(int i=0;i<num_vertices;i++){
		mbr->low[0] = min(mbr->low[0], p[i].x);
		mbr->high[0] = max(mbr->high[0], p[i].x);
		mbr->low[1] = min(mbr->low[1], p[i].y);
		mbr->high[1] = max(mbr->high[1], p[i].y);
	}
	return mbr;
}

void VertexSequence::fix(){
	if(num_vertices<=3){
		return;
	}
	int cur = 0;
	int next = 1;
	unordered_map<double, double> exist;
	exist[p[cur].x+p[cur].y] = p[cur].x;
	while(next<num_vertices){
		// next vertex appeared before, skip this one
		if(exist.find(p[next].x+p[next].y)!=exist.end() && exist[p[next].x+p[next].y] == p[next].x){
			next++;
			continue;
		}
		// not cllinear and shown first time, move cur forward
		if(!collinear(p[cur],p[next],p[(next+1)%num_vertices])){
			p[++cur] = p[next];
			// register the new one
			exist[p[cur].x+p[cur].y] = p[cur].x;
		}
		// the next vertex is valid
		next++;
	}
	exist.clear();
	num_vertices = cur+1;
}

size_t VertexSequence::encode(char *dest){
	assert(num_vertices>0);
	size_t encoded = 0;
	((long *)dest)[0] = num_vertices;
	encoded += sizeof(size_t);
	memcpy(dest+encoded,(char *)p,num_vertices*sizeof(Point));
	encoded += num_vertices*sizeof(Point);
	return encoded;
}

size_t VertexSequence::decode(char *source){
	size_t decoded = 0;
	num_vertices = ((size_t *)source)[0];
	assert(num_vertices>0);
	p = new Point[num_vertices];
	decoded += sizeof(size_t);
	memcpy((char *)p,source+decoded,num_vertices*sizeof(Point));
	decoded += num_vertices*sizeof(Point);
	return decoded;
}

size_t VertexSequence::get_data_size(){
	return sizeof(size_t)+num_vertices*sizeof(Point);
}


VertexSequence *VertexSequence::clone(){
	VertexSequence *ret = new VertexSequence(num_vertices);
	memcpy((void *)ret->p,(void *)p,sizeof(Point)*num_vertices);
	return ret;
}

/*
 *
 * containment test
 *
 * */
bool VertexSequence::contain(Point &point) {
	bool ret = false;
	for(int i = 0,j=num_vertices-1; i < num_vertices; j=i++) {
		if( ( (p[i].y >= point.y ) != (p[j].y >= point.y) ) &&
				(point.x <= (p[j].x - p[i].x) * (point.y - p[i].y) / (p[j].y - p[i].y) + p[i].x)){
			ret = !ret;
		}
	}
	return ret;
}

double VertexSequence::distance(Point &point, bool geography) {
	return point_to_segment_sequence_distance(point, p, num_vertices, geography);
}
