#include "../include/Ideal.h"

#include <iostream>
#include <iomanip>

int main(int argc, char** argv) {
    query_context global_ctx;
	vector<Ideal* > source = load_binary_file("/home/qmh/data/has_child.idl", global_ctx);
	RTree<Ideal *, double, 2, double> ideal_RTree;
	for(int i = 0; i < 30; i ++){
		auto p = source[i];
		ideal_RTree.Insert(p->getMBB()->low, p->getMBB()->high, p);
	}
	ideal_RTree.PrintTree();
	// long long int cnt = 0;
	// for(auto ideal : source){
	// 	ideal->rasterization(10);
	// 	for(int i = 0; i <= ideal->get_dimx(); i ++){
	// 		for(int j = 0; j <= ideal->get_dimy(); j ++){
	// 			int cur = ideal->get_id(i, j);
	// 			cout << ideal->get_num_sequences(cur) << " ";
	// 			for(int i = 0; i < ideal->get_num_sequences(cur); i ++){
	// 				auto rg = ideal->get_edge_sequence(ideal->get_offset(cur) + i);
	// 				cout << rg.first << " " << rg.second << "    ";
	// 			}
	// 			cout << endl;
	// 		}
	// 	}
	// 	cout << endl;
		// vector<int> boundary = ideal->get_pixels(BORDER);
		// for(auto cur : boundary){
		// 	cnt ++;
		// 	cout << ideal->get_num_sequences(cur) << " ";
		// 	for(int i = 0; i < ideal->get_num_sequences(cur); i ++){
		// 		auto rg = ideal->get_edge_sequence(ideal->get_offset(cur) + i);
		// 		cout << rg.second << " ";
		// 	}
		// 	cout << endl;
		// }
		// ideal->MyPolygon::print();
		// ideal->MyRaster::print();
	// }

	// cout << cnt << endl;



    return 0;
}





