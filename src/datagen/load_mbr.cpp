/*
 * Parser.cpp
 *
 *  Created on: May 9, 2022
 *      Author: teng
 */

#include "MyPolygon.h"
#include <boost/sort/sort.hpp>

using namespace std;

size_t load_binary_file(const char *path, box **mbrs, size_t **edge_num){
	if(!file_exist(path)){
		log("%s does not exist",path);
		return 0;
	}
	struct timeval start = get_cur_time();
	ifstream infile;
	infile.open(path, ios::in | ios::binary);
	size_t off;
	size_t num_polygons = 0;
	infile.seekg(0, infile.end);
	//seek to the first polygon
	infile.seekg(-sizeof(size_t), infile.end);
	infile.read((char *)&num_polygons, sizeof(size_t));
	assert(num_polygons>0 && "the file should contain at least one polygon");
	size_t *offsets = new size_t[num_polygons];

	infile.seekg(-sizeof(size_t)*(num_polygons+1), infile.end);
	infile.read((char *)offsets, sizeof(size_t)*num_polygons);

	log("loading %ld polygon from %s",num_polygons,path);

	size_t num_edges = 0;
	size_t data_size = 0;

	*edge_num = new size_t[num_polygons];
	*mbrs = new box[num_polygons];
	struct timeval timer = get_cur_time();

	for(size_t i=0;i<num_polygons;i++){
		if(get_time_elapsed(timer, false)>=200){
			log_refresh("loaded %f%%",i*100.0/num_polygons);
			timer = get_cur_time();
		}
		infile.seekg(offsets[i], infile.beg);
		MyPolygon *poly = MyPolygon::read_polygon_binary_file(infile);
		assert(poly);
		num_edges += poly->get_num_vertices();
		data_size += poly->get_data_size();
		(*mbrs)[i].low[0] = poly->getMBB()->low[0];
		(*mbrs)[i].high[0] = poly->getMBB()->high[0];
		(*mbrs)[i].low[1] = poly->getMBB()->low[1];
		(*mbrs)[i].high[1] = poly->getMBB()->high[1];
		(*edge_num)[i] = poly->get_num_vertices();
		delete poly;
	}
	infile.close();
	delete []offsets;
	logt("loaded %ld polygons each with %ld edges %ld MB", start, num_polygons,num_edges/num_polygons,data_size/1024/1024);
	return num_polygons;
}


int main(int argc, char** argv) {
	query_context global_ctx;
	global_ctx = get_parameters(argc, argv);

	box *boxes;
	size_t *edge_num;
	size_t box_num = load_binary_file(global_ctx.source_path.c_str(), &boxes, &edge_num);

	size_t *edge_num_tmp = new size_t[box_num];
	memcpy(edge_num_tmp, edge_num, box_num*sizeof(edge_num[0]));

	boost::sort::block_indirect_sort(edge_num_tmp, edge_num_tmp+box_num, std::greater<size_t>());

	size_t threshold = edge_num_tmp[box_num/4];
	size_t iter = 0;
	for(size_t i=0;i<box_num;i++){
		if(edge_num[i]<=threshold){
			boxes[iter++] = boxes[i];
		}
	}
	log("%ld smaller than %ld ",iter, threshold);
	write_file(global_ctx.target_path.c_str(), (char *)boxes, sizeof(box)*(box_num-box_num/4));

	delete []boxes;
	delete []edge_num;
	return 0;
}



