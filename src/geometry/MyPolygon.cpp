#include "../include/MyPolygon.h"

box *MyPolygon::getMBB(){
	if(mbr){
		return mbr;
	}
	mbr = boundary->getMBR();
	return mbr;
}

MyPolygon *MyPolygon::clone(){
	MyPolygon *polygon = new MyPolygon();
	polygon->boundary = boundary->clone();
	for(VertexSequence *t:holes){
		polygon->holes.push_back(t->clone());
	}
	polygon->getMBB();
	return polygon;
}

size_t MyPolygon::get_data_size(){
	size_t ds = 0;
	ds += sizeof(size_t);
	// for boundary
	ds += boundary->get_data_size();
	for(VertexSequence *vs:holes){
		ds += vs->get_data_size();
	}
	return ds;
}

VertexSequence *MyPolygon::read_vertices(const char *wkt, size_t &offset, bool clockwise){
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
	VertexSequence *vs = new VertexSequence(num_vertices);

	// read x/y
	for(int i=0;i<num_vertices;i++){
		vs->p[i].x = read_double(wkt, offset);
		vs->p[i].y = read_double(wkt, offset);
	}
	if(clockwise){
		if(!vs->clockwise()){
			vs->reverse();
		}
	}else{
		if(vs->clockwise()){
			vs->reverse();
		}
	}

	// read until the right parenthesis
	skip_space(wkt, offset);
	assert(wkt[offset++]==')');
	return vs;
}

MyPolygon *MyPolygon::read_polygon(const char *wkt, size_t &offset){
	MyPolygon *polygon = new MyPolygon();
	skip_space(wkt, offset);
	// left parentheses for the entire polygon
	assert(wkt[offset++]=='(');

	// read the vertices of the boundary polygon
	// the vertex must rotation in clockwise
	polygon->boundary = read_vertices(wkt, offset,false);
	if(polygon->boundary->clockwise()){
		polygon->boundary->reverse();
	}
	polygon->boundary->fix();
	skip_space(wkt, offset);
	//polygons as the holes of the boundary polygon
	while(wkt[offset]==','){
		offset++;
		VertexSequence *vc = read_vertices(wkt, offset,true);
		if(!vc->clockwise()){
			vc->reverse();
		}
		vc->fix();
		polygon->holes.push_back(vc);

		skip_space(wkt, offset);
	}
	assert(wkt[offset++]==')');
	polygon->getMBB();
	return polygon;
}

bool MyPolygon::validate_vertices(const char *wkt, size_t &offset, size_t &len){
	// read until the left parenthesis
	skip_space(wkt, offset);
	if(wkt[offset++]!='('){
		return false;
	}

	// count the number of vertices
	int cur_offset = offset;
	int num_vertices = 0;
	while(wkt[cur_offset++]!=')'){
		if(wkt[cur_offset]==','){
			num_vertices++;
		}
	}
	num_vertices++;

	// read x/y
	for(int i=0;i<num_vertices;i++){
		read_double(wkt, offset);
		read_double(wkt, offset);
	}

	// read until the right parenthesis
	skip_space(wkt, offset);
	if(wkt[offset++]!=')'){
		return false;
	}
	return true;
}

bool MyPolygon::validate_polygon(const char *wkt, size_t &offset, size_t &len){
	skip_space(wkt, offset);
	// left parentheses for the entire polygon
	if(wkt[offset++]!='('){
		return false;
	}

	// read the vertices of the boundary polygon
	// the vertex must rotation in clockwise
	if(!validate_vertices(wkt, offset, len)){
		return false;
	}
	skip_space(wkt, offset);
	//polygons as the holes of the boundary polygon
	while(wkt[offset]==','){
		offset++;
		if(!validate_vertices(wkt, offset, len)){
			return false;
		}
		skip_space(wkt, offset);
	}
	if(wkt[offset++]!=')'){
		return false;
	}
	return true;
}

MyPolygon *MyPolygon::gen_box(double min_x,double min_y,double max_x,double max_y){
	MyPolygon *mbr = new MyPolygon();
	mbr->boundary = new VertexSequence(5);
	mbr->boundary->p[0].x = min_x;
	mbr->boundary->p[0].y = min_y;
	mbr->boundary->p[1].x = max_x;
	mbr->boundary->p[1].y = min_y;
	mbr->boundary->p[2].x = max_x;
	mbr->boundary->p[2].y = max_y;
	mbr->boundary->p[3].x = min_x;
	mbr->boundary->p[3].y = max_y;
	mbr->boundary->p[4].x = min_x;
	mbr->boundary->p[4].y = min_y;
	return mbr;
}

MyPolygon *MyPolygon::gen_box(box &pix){
	return gen_box(pix.low[0],pix.low[1],pix.high[0],pix.high[1]);
}

void MyPolygon::print_without_head(bool print_hole, bool complete_ring){
	assert(boundary);
	cout<<"(";

	boundary->print(complete_ring);
	if(print_hole){
		for(VertexSequence *vs:this->holes){
			cout<<", ";
			vs->print(complete_ring);
		}
	}
	cout<<")";
}

void MyPolygon::print(bool print_id, bool print_hole){
	if(print_id){
		cout<<"id:\t"<<this->id<<endl;
	}
	cout<<"POLYGON";
	print_without_head(print_hole);
	cout<<endl;
}

PolygonMeta MyPolygon::get_meta(){
	PolygonMeta pmeta;
	pmeta.size = get_data_size();
	pmeta.num_vertices = get_num_vertices();
	pmeta.mbr = *getMBB();
	return pmeta;
}

size_t MyPolygon::encode(char *target){
	size_t encoded = 0;
	((size_t *)target)[0] = holes.size();
	encoded += sizeof(size_t); //saved one size_t for number of holes
	encoded += boundary->encode(target+encoded);
	for(VertexSequence *vs:holes){
		encoded += vs->encode(target+encoded);
	}
	return encoded;
}

size_t MyPolygon::decode(char *source){
	size_t decoded = 0;
	assert(!boundary);
	boundary = new VertexSequence();
	size_t num_holes = ((size_t *)source)[0];
	decoded += sizeof(size_t);
	decoded += boundary->decode(source+decoded);
	for(size_t i=0;i<num_holes;i++){
		VertexSequence *vs = new VertexSequence();
		decoded += vs->decode(source+decoded);
	}
	return decoded;
}

