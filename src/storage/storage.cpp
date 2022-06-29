/*
 * storage.cpp
 *
 *  Created on: Jun 27, 2022
 *      Author: teng
 */


#include "MyPolygon.h"


/*
 * in this file we define the .idl file format
 *
 * */

void dump_polygons_to_file(vector<MyPolygon *> polygons, const char *path){
	ofstream os;
	os.open(path, ios::out | ios::binary |ios::trunc);
	assert(os.is_open());

	size_t buffer_size = 100*1024*1024;
	char *data_buffer = new char[buffer_size];
	size_t data_size = 0;
	size_t curoffset = 0;
	vector<PolygonMeta> pmeta;
	pmeta.reserve(polygons.size());
	for(int i=0;i<polygons.size();i++){
		MyPolygon *p = polygons[i];
		if(p->get_data_size()+data_size>buffer_size){
			os.write(data_buffer, data_size);
			data_size = 0;
		}
		pmeta[i].offset = curoffset;
		pmeta[i].size = p->get_data_size();
		pmeta[i].mbr = *p->getMBB();
		pmeta[i].num_vertices = p->get_num_vertices();
		data_size += p->encode_to(data_buffer+data_size);
		curoffset += p->get_data_size();
	}

	// dump the rest polygon data
	if(data_size!=0){
		os.write(data_buffer, data_size);
	}
	// dump the meta data of the polygons
	os.write((char *)&pmeta[0], sizeof(PolygonMeta)*pmeta.size());
	size_t bs = polygons.size();
	os.write((char *)&bs, sizeof(size_t));
	os.close();

	pmeta.clear();
}

MyPolygon *read_polygon_binary_file(ifstream &infile){
	size_t num_holes = 0;
	infile.read((char *)&num_holes,sizeof(size_t));

	size_t num_vertices = 0;
	infile.read((char *)&num_vertices,sizeof(size_t));
	if(num_vertices==0){
		return NULL;
	}
	MyPolygon *poly = new MyPolygon();
	poly->boundary = new VertexSequence(num_vertices);
	infile.read((char *)poly->boundary->p,num_vertices*sizeof(Point));
	if(poly->boundary->clockwise()){
		poly->boundary->reverse();
	}

	for(int i=0;i<num_holes;i++){
		infile.read((char *)&num_vertices,sizeof(long));
		assert(num_vertices);
		VertexSequence *vs = new VertexSequence(num_vertices);
		infile.read((char *)vs->p,num_vertices*sizeof(Point));
		if(!vs->clockwise()){
			vs->reverse();
		}
		vs->fix();
		poly->holes.push_back(vs);
	}
	return poly;
}


// idx starting from 0
MyPolygon *load_binary_file_single(const char *path, query_context ctx, int idx){
	ifstream infile;
	infile.open(path, ios::in | ios::binary);

	size_t num_polygons_infile;
	infile.seekg(-sizeof(size_t), infile.end);
	infile.read((char *)&num_polygons_infile, sizeof(size_t));
	assert(idx<num_polygons_infile && "the idx must smaller than the polygon number ");

	PolygonMeta pmeta;
	infile.seekg(-sizeof(size_t) - sizeof(PolygonMeta)*(num_polygons_infile-idx), infile.end);
	infile.read((char *)&pmeta, sizeof(PolygonMeta));

	char *buffer = new char[pmeta.size];

	infile.seekg(pmeta.offset, infile.beg);
	infile.read(buffer, pmeta.size);

	MyPolygon *poly = new MyPolygon();
	poly->decode_from(buffer);

	delete []buffer;
	infile.close();
	return poly;
}


typedef struct{
	ifstream *infile;
	size_t offset;
	size_t poly_size;
	size_t load(char *buffer){
		infile->seekg(offset, infile->beg);
		infile->read(buffer, poly_size);
		return poly_size;
	}
}load_holder;


const size_t buffer_size = 10*1024*1024;

void *load_unit(void *arg){
	query_context *ctx = (query_context *)arg;
	vector<load_holder *> *jobs = (vector<load_holder *> *)ctx->target;
	vector<MyPolygon *> *global_polygons = (vector<MyPolygon *> *)ctx->target2;

	char *buffer = new char[buffer_size];
	vector<MyPolygon *> polygons;
	while(ctx->next_batch(1)){
		for(int i=ctx->index;i<ctx->index_end;i++){
			load_holder *lh = (*jobs)[i];
			ctx->global_ctx->lock();
			size_t poly_size = lh->load(buffer);
			ctx->global_ctx->unlock();
			size_t off = 0;
			while(off<poly_size){
				MyPolygon *poly = new MyPolygon();
				off += poly->decode_from(buffer+off);
				polygons.push_back(poly);
			}
			ctx->report_progress(1);
		}
	}
	delete buffer;
	ctx->global_ctx->lock();
	global_polygons->insert(global_polygons->end(), polygons.begin(), polygons.end());
	ctx->global_ctx->unlock();
	polygons.clear();
	return NULL;
}

vector<MyPolygon *> load_binary_file(const char *path, query_context &ctx, bool sample){
	vector<MyPolygon *> polygons;
	if(!file_exist(path)){
		log("%s does not exist",path);
		exit(0);
	}
	struct timeval start = get_cur_time();

	ifstream infile;
	infile.open(path, ios::in | ios::binary);
	size_t num_polygons_infile = 0;
	infile.seekg(0, infile.end);
	//seek to the first polygon
	infile.seekg(-sizeof(size_t), infile.end);
	infile.read((char *)&num_polygons_infile, sizeof(size_t));
	assert(num_polygons_infile>0 && "the file should contain at least one polygon");

	PolygonMeta *pmeta = new PolygonMeta[num_polygons_infile];
	infile.seekg(-sizeof(size_t)-sizeof(PolygonMeta)*num_polygons_infile, infile.end);
	infile.read((char *)pmeta, sizeof(PolygonMeta)*num_polygons_infile);
	// the last one is the end
	size_t num_polygons = min(num_polygons_infile, ctx.max_num_polygons);

	logt("loading %ld polygon from %s",start, num_polygons,path);
	// organizing tasks
	vector<load_holder *> tasks;
	size_t cur = 0;
	while(cur<num_polygons){
		size_t end = cur+1;

		while(end<num_polygons &&
				pmeta[end].offset - pmeta[cur].offset + pmeta[end].size < buffer_size){
			end++;
		}
		load_holder *lh = new load_holder();
		lh->infile = &infile;
		lh->offset = pmeta[cur].offset;
		if(end<num_polygons){
			lh->poly_size = pmeta[end].offset - pmeta[cur].offset;
		}else{
			lh->poly_size = pmeta[end-1].offset - pmeta[cur].offset + pmeta[end-1].size;
		}
		tasks.push_back(lh);
		cur = end;
	}

	logt("packed %ld tasks", start, tasks.size());

	query_context global_ctx;
	global_ctx.target_num = tasks.size();
	pthread_t threads[global_ctx.num_threads];
	query_context myctx[global_ctx.num_threads];
	for(int i=0;i<global_ctx.num_threads;i++){
		myctx[i] = global_ctx;
		myctx[i].thread_id = i;
		myctx[i].global_ctx = &global_ctx;
		myctx[i].target = (void *)&tasks;
		myctx[i].target2 = (void *)&polygons;
	}
	for(int i=0;i<global_ctx.num_threads;i++){
		pthread_create(&threads[i], NULL, load_unit, (void *)&myctx[i]);
	}

	for(int i = 0; i < global_ctx.num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}
	infile.close();
	delete []pmeta;
	for(load_holder *lh:tasks){
		delete lh;
	}
	logt("loaded %ld polygons", start, polygons.size());
	return polygons;
}

size_t load_polygonmeta_from_file(const char *path, PolygonMeta **pmeta){
	ifstream infile;
	infile.open(path, ios::in | ios::binary);
	size_t num_polygons_infile = 0;
	infile.seekg(0, infile.end);
	//seek to the first polygon
	infile.seekg(-sizeof(size_t), infile.end);
	infile.read((char *)&num_polygons_infile, sizeof(size_t));
	assert(num_polygons_infile>0 && "the file should contain at least one polygon");

	*pmeta = new PolygonMeta[num_polygons_infile];
	infile.seekg(-sizeof(size_t)-sizeof(PolygonMeta)*num_polygons_infile, infile.end);
	infile.read((char *)*pmeta, sizeof(PolygonMeta)*num_polygons_infile);

	return num_polygons_infile;
}

size_t load_mbr_from_file(const char *path, box **mbrs){

	if(!file_exist(path)){
		log("%s is empty",path);
		exit(0);
	}

	PolygonMeta *pmeta;
	size_t num_polygons = load_polygonmeta_from_file(path, &pmeta);
	*mbrs = new box[num_polygons];
	for(size_t i=0;i<num_polygons;i++){
		(*mbrs)[i] = pmeta[i].mbr;
	}
	delete []pmeta;
	return num_polygons;
}

size_t load_points_from_path(const char *path, Point **points){
	size_t fsize = file_size(path);
	if(fsize<=0){
		log("%s is empty",path);
		exit(0);
	}
	size_t target_num = fsize/sizeof(Point);
	log_refresh("start loading %ld points",target_num);

	*points = new Point[target_num];
	ifstream infile(path, ios::in | ios::binary);
	infile.read((char *)*points, fsize);
	infile.close();
	return target_num;
}


