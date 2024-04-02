#include "../include/MyPolygon.h"
#include "../include/MyRaster.h"

box *MyPolygon::getMBB(){
	if(mbr){
		return mbr;
	}
	mbr = boundary->getMBR();
	return mbr;
}

box *MyPolygon::getMER(query_context *ctx){

	if(mer){
		return mer;
	}

	MyRaster *ras = new MyRaster();
	ras->rasterization(boundary, ctx->vpr);
	vector<int> interiors = ras->get_pixels(IN);
	if(interiors.size()==0){
		return NULL;
	}
	int loops = ctx->mer_sample_round;
	int maxcount = 0;
	box *max_mer = NULL;
	while(loops-->0){
		int sample = get_rand_number(interiors.size())-1;
		box *curmer = ras->extractMER(interiors[sample]);
		if(max_mer){
			if(max_mer->area()<curmer->area()){
				delete max_mer;
				max_mer = curmer;
			}else{
				delete curmer;
			}
		}else{
			max_mer = curmer;
		}
	}
	interiors.clear();
	mer = max_mer;
	delete ras;
	return mer;
}

VertexSequence *MyPolygon::get_convex_hull(){
	if(convex_hull==NULL){
		convex_hull = boundary->convexHull();
	}
	return convex_hull;
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

size_t MyPolygon::get_rtree_size(){
	if(rtree){
		int count = rtree->node_count();
		return count*4*8;
	}else{
		return 0;
	}
}

void MyPolygon::triangulate(){
	assert(!boundary->clockwise());
	if(triangle_num > 0){
		return;
	}
	vector<Vertex *> polyline = boundary->pack_to_polyline();
	assert(polyline.size() > 0);
	CDT *cdt = new CDT(polyline);
	cdt->Triangulate();
	vector<Triangle *> tri = cdt->GetTriangles();
	triangles = new Point[3*tri.size()];
	for(int i=0;i<tri.size();i++){
		triangles[i*3].x = tri[i]->point(0)->x;
		triangles[i*3].y = tri[i]->point(0)->y;
		triangles[i*3+1].x = tri[i]->point(1)->x;
		triangles[i*3+1].y = tri[i]->point(1)->y;
		triangles[i*3+2].x = tri[i]->point(2)->x;
		triangles[i*3+2].y = tri[i]->point(2)->y;
	}
	triangle_num = tri.size();
	delete cdt;
	for(Vertex *v:polyline){
		delete v;
	}
	polyline.clear();
}

void MyPolygon::build_rtree(){
	triangulate();

	assert(triangles && triangle_num>0);
	if(rtree){
		return;
	}

	RTree<Point *, double, 2, double> *rtree_tmp = new RTree<Point *, double, 2, double>();

	for(int i=0;i<triangle_num;i++){
		box pix;
		Point *ps = triangles+3*i;
		for(int i=0;i<3;i++){
			pix.update(ps[i]);
		}
		rtree_tmp->Insert(pix.low, pix.high, ps);
	}
	mbr = getMBB();
	rtree = new RTNode();
	rtree->low[0] = mbr->low[0];
	rtree->low[1] = mbr->low[1];
	rtree->high[0] = mbr->high[0];
	rtree->high[1] = mbr->high[1];

	rtree_tmp->construct_pixel(rtree);
	assert(rtree->validate());
	delete rtree_tmp;
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

bool MyPolygon::contain(Point &p){
    return boundary->contain(p);
}

bool MyPolygon::contain_rtree(RTNode *node, Point &p, query_context *ctx){
	if(!node->contain(p)){
		return false;
	}
	if(node->is_leaf()){
		assert(node->node_element);
		Point *triangle = (Point *)node->node_element;
		struct timeval start = get_cur_time();
		bool ret = false;
		for(int i=0;i<=2;i++){
			Point *start = triangle+i;
			Point *end = triangle+(i+1)%3;
			if((start->y >= p.y ) != (end->y >= p.y)){
				double xint = (end->x - start->x) * (p.y - start->y)/ (end->y - start->y) + start->x;
				if(p.x <= xint){
					ret = !ret;
				}
			}
		}
		return !ret;
	}
	for(RTNode *ch:node->children){
		if(contain_rtree(ch,p,ctx)){
			return true;
		}
	}
	return false;
}

bool MyPolygon::contain(Point &p, query_context *ctx, bool profile){
	// check the maximum enclosed rectangle (MER)
	if(mer&&mer->contain(p)){
		return true;
	}
	// check the convex hull
	if(convex_hull&&!convex_hull->contain(p)){
		return false;
	}

	// refinement step
	timeval start = get_cur_time();
	bool contained = false;

	// todo for test only, remove in released version
	if(ctx->perform_refine)
	{
		// if(rtree){
		// 	contained = contain_rtree(rtree,p,ctx);
		// }else{
			contained = contain(p);
		// }
	}
	return contained;
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

