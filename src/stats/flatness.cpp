/*
 * partstat.cpp
 *
 *  Created on: Jun 4, 2020
 *      Author: teng
 */
/*
 * Parser.cpp
 *
 *  Created on: May 9, 2020
 *      Author: teng
 */

#include "../geometry/MyPolygon.h"
#include <fstream>
#include "../index/RTree.h"
#include <queue>
#include <boost/program_options.hpp>
#include "cuda/mygpu.h"

namespace po = boost::program_options;
using namespace std;

int element_size = 100;

// some shared parameters
pthread_mutex_t poly_lock;
bool stop = false;

pthread_mutex_t report_lock;
long query_count = 0;


queue<queue_element *> shared_queue;


class partition_context{
public:
	int vpr = 10;
	vector<double> flatneses;
	vector<int> borders;
	vector<int> ins;
	vector<int> outs;
	vector<int> crosses_count;
	size_t total_num_vertices = 0;
};




void *partition(void *args){
	partition_context *ctx = (partition_context *)args;
	vector<queue_element *> elems;
	int totalflat = ctx->flatneses.size();
	int border[totalflat];
	int in[totalflat];
	int out[totalflat];
	int cross_count[totalflat];
	size_t num_vertices = 0;
	for(int i=0;i<totalflat;i++){
		border[i] = 0;
		in[i] = 0;
		out[i] = 0;
		cross_count[i] = 0;
	}
	while(!stop||shared_queue.size()>0){
		queue_element *elem = NULL;
		pthread_mutex_lock(&poly_lock);
		if(!shared_queue.empty()){
			elem = shared_queue.front();
			shared_queue.pop();
		}
		pthread_mutex_unlock(&poly_lock);
		// queue is empty but not stopped
		if(!elem){
			usleep(20);
			continue;
		}
		for(MyPolygon *poly:elem->polys){
			num_vertices += poly->get_num_vertices();
			for(int k=0;k<totalflat;k++){
				vector<vector<Pixel>> partitions = poly->partition(ctx->vpr, ctx->flatneses[k]);
				for(int i=0;i<partitions.size();i++){
					for(int j=0;j<partitions[0].size();j++){
						if(partitions[i][j].status==BORDER){
							border[k]++;
						}else if(partitions[i][j].status==IN){
							in[k]++;
						}else if(partitions[i][j].status==OUT){
							out[k]++;
						}
						cross_count[k]+=partitions[i][j].crosses.size();
					}
				}
				poly->reset_grid_partition();
			}
		}
		delete elem;
	}

	pthread_mutex_lock(&report_lock);
	for(int i=0;i<totalflat;i++){
		ctx->borders[i]+=border[i];
		ctx->ins[i]+=in[i];
		ctx->outs[i]+=out[i];
		ctx->crosses_count[i] += cross_count[i];
	}
	ctx->total_num_vertices += num_vertices;

	pthread_mutex_unlock(&report_lock);

	return NULL;
}



int main(int argc, char** argv) {
	string source_path;
	int num_threads = get_num_threads();
	partition_context ctx;

	double max_flatness = 10;
	double flatness_step = 0.5;
	po::options_description desc("query usage");
	desc.add_options()
		("help,h", "produce help message")
		("vpr,v", po::value<int>(&ctx.vpr), "number of vertices per raster")
		("source_path,s", po::value<string>(&source_path), "the path to the source")
		("max_flatness", po::value<double>(&max_flatness), "the maximum flatness value")
		("flatness_step", po::value<double>(&flatness_step), "the flatness step")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	if (vm.count("help")) {
		cout << desc << "\n";
		return 0;
	}
	po::notify(vm);

	for(double s=max_flatness;s>1;s-=flatness_step){
		ctx.flatneses.push_back(1/s);
	}
	ctx.flatneses.push_back(1);
	for(double s=1+flatness_step;s<=max_flatness;s+=flatness_step){
		ctx.flatneses.push_back(s);
	}

	for(int i=0;i<ctx.flatneses.size();i++){
		ctx.outs.push_back(0);
		ctx.ins.push_back(0);
		ctx.borders.push_back(0);
		ctx.crosses_count.push_back(0);
	}


	timeval start = get_cur_time();

	pthread_t threads[num_threads];

	for(int i=0;i<num_threads;i++){
		pthread_create(&threads[i], NULL, partition, (void *)&ctx);
	}

	ifstream is;
	is.open(source_path.c_str(), ios::in | ios::binary);
	int eid = 0;
	queue_element *elem = new queue_element(eid++);
	int loaded = 0;
	size_t pid = 0;
	while(!is.eof()){
		MyPolygon *poly = MyPolygon::read_polygon_binary_file(is);
		if(!poly){
			continue;
		}
		poly->setid(pid++);
		elem->polys.push_back(poly);
		if(elem->polys.size()==element_size){
			while(shared_queue.size()>10*num_threads){
				usleep(20);
			}
			pthread_mutex_lock(&poly_lock);
			shared_queue.push(elem);
			pthread_mutex_unlock(&poly_lock);
			elem = new queue_element(eid++);
		}
		if(++loaded%100000==0){
			log("partitioned %d polygons",loaded);
		}
	}
	if(elem->polys.size()>0){
		pthread_mutex_lock(&poly_lock);
		shared_queue.push(elem);
		pthread_mutex_unlock(&poly_lock);
	}else{
		delete elem;
	}
	stop = true;

	for(int i = 0; i < num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}
	logt("partitioned %d polygons",start,loaded);
	for(int i=0;i<ctx.flatneses.size();i++){
		int total = ctx.ins[i]+ctx.outs[i]+ctx.borders[i];
//		if(total==0){
//			cout<<min_flatness+i*flatness_step<<endl;
//		}
		assert(total>0);

		printf("%f\t%f\t%f\n",ctx.flatneses[i],1.0*ctx.borders[i]/total,0.5*ctx.crosses_count[i]/ctx.total_num_vertices);
	}


	return 0;
}



