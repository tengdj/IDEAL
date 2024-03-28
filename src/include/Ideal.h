#ifndef IDEAL_H
#define IDEAL_H

#include "MyPolygon.h"
#include "MyRaster.h"

enum Direction{
	HORIZONTAL = 0,
	VERTICAL = 1
};

enum cross_type{
	ENTER = 0,
	LEAVE = 1
};

class cross_info{
public:
	cross_type type;
	int edge_id;
	cross_info(cross_type t, int e){
		type = t;
		edge_id = e;
	}
};



class Grid_line{
	uint16_t *offset = nullptr;
	double *intersection_nodes = nullptr;

	size_t num_grid_lines = 0;
	size_t num_crosses = 0;
public:
	Grid_line() = default;
	Grid_line(int size);
	~Grid_line();
	void init_intersection_node(int num_nodes);
	int get_num_nodes(int y) {return offset[y + 1] - offset[y];}
	void add_node(int idx, double x) {intersection_nodes[idx] = x;}

	void set_num_crosses(size_t x) {num_crosses = x;}
	size_t get_num_crosses() {return num_crosses;}
	void set_offset(int id, int idx) {offset[id] = idx;}
	uint16_t get_offset(int id) {return offset[id];}
	double get_intersection_nodes(int id) {return intersection_nodes[id];}
};


class Ideal : public MyPolygon, public MyRaster{
	uint16_t *offset = nullptr;
	pair<uint16_t, uint16_t> *edge_sequences = nullptr;   
	Grid_line *horizontal;
	Grid_line *vertical;

	int len_edge_sequences = 0;

    pthread_mutex_t ideal_partition_lock;
	void init_pixels();
	void evaluate_edges();
	void scanline_reandering();

public:
    Ideal(){
        pthread_mutex_init(&ideal_partition_lock, NULL);
    }
    void rasterization(int vertex_per_raster);
	void rasterization();

	void set_offset(int id, int idx){offset[id] = idx;}
	uint16_t get_offset(int id) {return offset[id];}
	void process_pixels_null(int x, int y);
	void init_edge_sequences(int num_edge_seqs);
	void add_edge(int idx, int start, int end);
	pair<uint16_t, uint16_t> get_edge_sequence(int idx){return edge_sequences[idx];}
	uint16_t get_num_sequences(int id);
	// double get_possible_min(Point &p, int center, int step, bool geography = true);
	void process_crosses(map<int, vector<cross_info>> edge_info);
	void process_intersection(map<int, vector<double>> edge_intersection, Direction direction);
	int count_intersection_nodes(Point &p);



	// statistic collection
	int get_num_border_edge();
	int num_edges_covered(int id);
	// size_t get_num_gridlines();
	size_t get_num_crosses();
	// double get_num_intersection();
	

	// query functions
	bool contain(Point &p, query_context *ctx, bool profile = false);
	bool contain(Ideal *target, query_context *ctx, bool profile = false);
	// bool intersect(MyPolygon *target, query_context *ctx);
	// double distance(Point &p, query_context *ctx, bool profile = true);
	// double distance(MyPolygon *target, query_context *ctx);
	// double distance(MyPolygon *target, int pix, query_context *ctx, bool profile = true);
};

//utility functions
void process_rasterization(query_context *ctx);
void preprocess(query_context *gctx);

// storage related functions
vector<Ideal *> load_binary_file(const char *path, query_context &ctx);

#endif // IDEAL_H