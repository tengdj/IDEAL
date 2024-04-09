#include "Ideal.h"
#include "query_context.h"

int main(int argc, char** argv) {
    query_context global_ctx;
	global_ctx.num_threads = 1;
	vector<Ideal* > source = load_binary_file("/home/qmh/data/has_child.idl", global_ctx);
	long long int cnt = 0;
	for(auto ideal : source){
		ideal->rasterization(10);
		for(int i = 0; i <= ideal->get_dimx(); i ++){
			for(int j = 0; j <= ideal->get_dimy(); j ++){
				int cur = ideal->get_id(i, j);
				cout << ideal->get_num_sequences(cur) << " ";
				for(int i = 0; i < ideal->get_num_sequences(cur); i ++){
					auto rg = ideal->get_edge_sequence(ideal->get_offset(cur) + i);
					cout << rg.first << " " << rg.second << "    ";
				}
				cout << endl;
			}
		}
		cout << endl;
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
	}

	cout << cnt << endl;



    return 0;
}
