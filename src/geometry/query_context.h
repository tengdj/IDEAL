/*
 * query_context.h
 *
 *  Created on: Jan 1, 2021
 *      Author: teng
 */

#ifndef SRC_GEOMETRY_QUERY_CONTEXT_H_
#define SRC_GEOMETRY_QUERY_CONTEXT_H_
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include "Point.h"
#include "Pixel.h"

using namespace std;
class MyPolygon;


enum QueryType{
    contain = 0,
    distance = 1,
    within = 2
};

class query_context{
public:
	//configuration
	int thread_id = 0;
	int num_threads = 0;
	int vpr = 10;
	int vpr_end = 10;
	bool use_grid = false;
	bool use_qtree = false;
	bool use_convex_hull = false;
	bool use_mer = false;
	bool use_triangulate = false;
	int mer_sample_round = 20;
	bool perform_refine = true;
	bool gpu = false;
	bool collect_latency = false;
	float sample_rate = 1.0;
	int small_threshold = 500;
	int big_threshold = 1000000;
	bool sort_polygons = false;
	int report_gap = 1000000;
	int distance_buffer_size = 10;
	QueryType query_type = QueryType::contain;
	string source_path;
	string target_path;

	string valid_path;

	//query target, for temporary use
	void *target = NULL;
	void *target2 = NULL;
	query_context *global_ctx = NULL;

	//shared staff
	size_t index = 0;
	size_t index_end = 0;
	pthread_mutex_t lock;


	//result
	double distance = 0;
	bool filter_checked_only = false;

	//query statistic
	size_t found = 0;
	size_t query_count = 0;
	size_t checked_count = 0;
	size_t refine_count = 0;

	size_t pixel_checked = 0;
	size_t border_checked = 0;
	size_t edges_checked = 0;
	double pixel_check_time = 0;
	double edges_check_time = 0;
	double check_time = 0;


	// temporary storage for query processing
	vector<MyPolygon *> source_polygons;
	vector<MyPolygon *> target_polygons;
	double *points = NULL;
	size_t target_num = 0;


	map<int, int> vertex_number;
	map<int, double> latency;


	// functions
	query_context();
	~query_context();
	query_context(query_context &t);
	query_context operator + (query_context const &obj);
	query_context& operator=(query_context const &obj);

	void report_latency(int num_v, double latency);
	void load_points();
	void report_progress();
	void merge_global();
	bool next_batch(int batch_num=1);
	void reset_stats(){
		//query statistic
		found = 0;
		query_count = 0;
		checked_count = 0;
		refine_count = 0;
		pixel_checked = 0;
		border_checked = 0;
		edges_checked = 0;
		pixel_check_time = 0;
		edges_check_time = 0;
		check_time = 0;
		index = 0;
	}

	bool is_within_query(){
		return query_type == QueryType::within;
	}

	bool is_contain_query(){
		return query_type == QueryType::contain;
	}
	void print_stats();
};


query_context get_parameters(int argc, char **argv);


#endif /* SRC_GEOMETRY_QUERY_CONTEXT_H_ */
