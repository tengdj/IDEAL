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
	int id = 0;
	int vpr = 10;
	double min_flatness = 0.1;
	double max_flatness = 2.0;
	double flatness_step = 0.1;
};

vector<int> borders;
vector<int> ins;
vector<int> outs;


void *partition(void *args){
	partition_context *ctx = (partition_context *)args;
	pthread_mutex_lock(&poly_lock);
	log("thread %d is started",ctx->id);
	pthread_mutex_unlock(&poly_lock);
	vector<queue_element *> elems;
	int totalflat = (ctx->max_flatness-ctx->min_flatness)/ctx->flatness_step;
	int border[totalflat];
	int in[totalflat];
	int out[totalflat];
	for(int i=0;i<totalflat;i++){
		border[i] = 0;
		in[i] = 0;
		out[i] = 0;
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
			for(int k=0;k<totalflat;k++){
				double flatness = ctx->min_flatness+k*ctx->flatness_step;
				vector<vector<Pixel>> partitions = poly->partition(ctx->vpr, flatness);
				for(int i=0;i<partitions.size();i++){
					for(int j=0;j<partitions[0].size();j++){
						MyPolygon *m = partitions[i][j].to_polygon();
						if(partitions[i][j].status==BORDER){
							border[k]++;
						}else if(partitions[i][j].status==IN){
							in[k]++;
						}else if(partitions[i][j].status==OUT){
							out[k]++;
						}
					}
				}
				poly->reset_partition();
			}
		}
		delete elem;
	}

	pthread_mutex_lock(&report_lock);
	for(int i=0;i<totalflat;i++){
		borders[i]+=border[i];
		ins[i]+=in[i];
		outs[i]+=out[i];
	}

	pthread_mutex_unlock(&report_lock);

	return NULL;
}



int main(int argc, char** argv) {
	string source_path;
	int num_threads = get_num_threads();
	int vpr = 10;
	double min_flatness = 0.1;
	double max_flatness = 2.0;
	double flatness_step = 0.1;
	po::options_description desc("query usage");
	desc.add_options()
		("help,h", "produce help message")
		("vpr,v", po::value<int>(&vpr), "number of vertices per raster")
		("source_path,s", po::value<string>(&source_path), "the path to the source")
		("min_flatness", po::value<double>(&min_flatness), "the minimum flatness value")
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
	assert(min_flatness>0&&max_flatness>min_flatness&&flatness_step>0);
	int totalflat = (max_flatness-min_flatness)/flatness_step;
	for(int i=0;i<totalflat;i++){
		outs.push_back(0);
		ins.push_back(0);
		borders.push_back(0);
	}


	timeval start = get_cur_time();

	pthread_t threads[num_threads];
	partition_context ctx[num_threads];
	for(int i=0;i<num_threads;i++){
		ctx[i].id = i;
		ctx[i].vpr = vpr;
		ctx[i].min_flatness = min_flatness;
		ctx[i].max_flatness = max_flatness;
		ctx[i].flatness_step = flatness_step;
	}


	for(int i=0;i<num_threads;i++){
		pthread_create(&threads[i], NULL, partition, (void *)&ctx[i]);
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
	for(int i=0;i<totalflat;i++){
		int total = ins[i]+outs[i]+borders[i];
		if(total==0){
			cout<<min_flatness+i*flatness_step<<endl;
		}
		assert(total>0);

		printf("%f\t%f\t%f\t%f\n",min_flatness+i*flatness_step,1.0*borders[i]/total,1.0*ins[i]/total,1.0*outs[i]/total);
	}


	return 0;
}



