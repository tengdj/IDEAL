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
	box *boxes = new box[polygons.size()];
	size_t *offsets = new size_t[polygons.size()];
	size_t curoffset = 0;
	for(int i=0;i<polygons.size();i++){
		MyPolygon *p = polygons[i];
		if(p->get_data_size()+data_size>buffer_size){
			os.write(data_buffer, data_size);
			data_size = 0;
		}
		data_size += p->encode_to(data_buffer+data_size);
		offsets[i] = curoffset;
		boxes[i] = *p->getMBB();
		curoffset += p->get_data_size();
	}

	// dump the rest polygon data
	if(data_size!=0){
		os.write(data_buffer, data_size);
	}
	// dump the bounding boxes
	os.write((char *)boxes, sizeof(box)*polygons.size());
	// dump the offsets of the polygons
	os.write((char *)offsets, sizeof(size_t)*polygons.size());
	size_t bs = polygons.size();
	os.write((char *)&bs, sizeof(size_t));
	os.close();

	delete []boxes;
	delete []offsets;
}

//MyPolygon *read_polygon_binary_file(ifstream &infile){
//	size_t num_holes = 0;
//	infile.read((char *)&num_holes,sizeof(size_t));
//
//	long num_vertices = 0;
//	infile.read((char *)&num_vertices,sizeof(long));
//	if(num_vertices==0){
//		return NULL;
//	}
//	MyPolygon *poly = new MyPolygon();
//	poly->boundary = new VertexSequence(num_vertices);
//	infile.read((char *)poly->boundary->p,num_vertices*sizeof(Point));
//	if(poly->boundary->clockwise()){
//		poly->boundary->reverse();
//	}
//
//	for(int i=0;i<num_holes;i++){
//		infile.read((char *)&num_vertices,sizeof(long));
//		assert(num_vertices);
//		VertexSequence *vs = new VertexSequence(num_vertices);
//		infile.read((char *)vs->p,num_vertices*sizeof(Point));
//		if(!vs->clockwise()){
//			vs->reverse();
//		}
//		vs->fix();
//		poly->holes.push_back(vs);
//	}
//	return poly;
//}


// idx starting from 0
MyPolygon *load_binary_file_single(const char *path, query_context ctx, int idx){
	ifstream infile;
	infile.open(path, ios::in | ios::binary);

	size_t num_polygons_infile;
	infile.seekg(-sizeof(size_t), infile.end);
	infile.read((char *)&num_polygons_infile, sizeof(size_t));
	assert(idx<num_polygons_infile && "the idx must smaller than the polygon number ");

	size_t offset;
	infile.seekg(-sizeof(size_t)*(num_polygons_infile+1-idx), infile.end);
	infile.read((char *)&offset, sizeof(size_t));

	size_t poly_size = 0;
	// the last one
	if(idx == num_polygons_infile-1){
		poly_size = file_size(path) - sizeof(size_t)*(num_polygons_infile+1) - offset;
	}else{
		size_t next_offset;
		infile.read((char *)&next_offset, sizeof(size_t));
		poly_size = next_offset - offset;
	}

	char *buffer = new char[poly_size];

	infile.seekg(offset,infile.beg);
	infile.read(buffer, poly_size);

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

	size_t *offsets = new size_t[num_polygons_infile+1];
	infile.seekg(-sizeof(size_t)*(num_polygons_infile+1), infile.end);
	infile.read((char *)offsets, sizeof(size_t)*num_polygons_infile);
	// the last one is the end
	offsets[num_polygons_infile] = file_size(path) - sizeof(size_t)*(num_polygons_infile+1);
	size_t num_polygons = min(num_polygons_infile, ctx.max_num_polygons);

	logt("loading %ld polygon from %s",start, num_polygons,path);

	// organizing tasks
	vector<load_holder *> tasks;
	size_t cur = 0;
	while(cur<num_polygons){
		size_t end = cur+1;
		while(end<=num_polygons){
			if(end==num_polygons ||
					offsets[end]-offsets[cur]>=buffer_size){
				load_holder *lh = new load_holder();
				lh->infile = &infile;
				lh->offset = offsets[cur];
				lh->poly_size = offsets[end - (end<num_polygons)] - offsets[cur];
				cur = end - (end<num_polygons);
				tasks.push_back(lh);
				break;
			}
			end++;
		}
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
	delete []offsets;
	for(load_holder *lh:tasks){
		delete lh;
	}
	logt("loaded %ld polygons", start, polygons.size());
	return polygons;
}

size_t load_boxes_from_file(const char *path, box **mbrs){
	size_t fsize = file_size(path);
	if(fsize<=0){
		log("%s is empty",path);
		exit(0);
	}
	size_t target_num = fsize/sizeof(box);
	log_refresh("start loading %ld MBRs",target_num);
	*mbrs = new box[target_num];

	ifstream infile(path, ios::in | ios::binary);
	infile.read((char *)*mbrs, fsize);
	infile.close();
	return target_num;
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


