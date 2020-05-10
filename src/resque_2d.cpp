
#include "resque_2d.hpp"
#include "MyPolygon.h"

/* 
 * RESQUE processing engine v3.0
 *   It supports spatial join and nearest neighbor query with different predicates
 *   1) parseParameters
 *   2) readCacheFile - metadata such as partition schemata
 *   3) for every input line in the current tile
 *         an input line represents an object
 *         save geometry and original data in memory
 *         execute join operation when finish reading a tile
 *   4) Join operation between 2 sets or a single set
 *         build Rtree index on the second data set
 *         for every object in the first data set
 *            using Rtree index of the second data set
 *              check for MBR/envelope intersection
 *              output the pair result or save pair statistics
 *   5) Output final statistics (opt)
 *   Requirement (input files): see the Wiki
 * */

using namespace geos;
using namespace geos::io;
using namespace geos::geom;
using namespace geos::operation::buffer;
using namespace geos::operation::distance;
using namespace std;
using namespace SpatialIndex;

#define QUEUE_SIZE 100
#define MAX_THREAD_NUM 100

// some shared parameters
string processing_line[MAX_THREAD_NUM];
bool is_working[MAX_THREAD_NUM];
pthread_mutex_t line_lock;
pthread_mutex_t output_lock;
bool stop = false;
int k = 100;

class geometry_wrapper{
public:
	double area;
	Geometry *geom=NULL;
	string raw;
	~geometry_wrapper(){
		if(geom!=NULL){
			delete geom;
		}
	}
};

vector<geometry_wrapper *> global_top;

void update(vector<geometry_wrapper *> &top, geometry_wrapper *geo){
	if(top.size()<k||top[k-1]->area<geo->area){
		if(top.size()==k){
			delete top[k-1];
			top.erase(top.begin()+k-1);
		}

		int insert_index = 0;
		for(insert_index=0;insert_index<top.size();insert_index++){
			if(top[insert_index]->area<geo->area){
				break;
			}
		}
		top.insert(top.begin()+insert_index, geo);
	}else{
		delete geo;
	}
}

MyPolygon *max_poly = NULL;


void *process_wkt(void *args){

	int id = *(int *)args;
	pthread_mutex_lock(&output_lock);
	log("thread %d is started", id);
	pthread_mutex_unlock(&output_lock);
	MyPolygon *local_max_poly = NULL;

	// Reading from the cache file
	/* Parsing polygon input */
	while (!stop||is_working[id]) {
		if(!is_working[id]){
			usleep(10);
			continue;
		}
		MyMultiPolygon *mp = new MyMultiPolygon(processing_line[id].c_str());
		vector<MyPolygon *> polygons = mp->get_polygons();
		for(MyPolygon *p:polygons){
			if(!local_max_poly){
				local_max_poly = p->clone();
			}else if(p->num_boundary_vertices()>local_max_poly->num_boundary_vertices()){
				//cout<<"err"<<endl;
				delete local_max_poly;
				local_max_poly = p->clone();
				//cout<<"hh"<<endl;
			}
		}

		delete mp;

		processing_line[id].clear();
		is_working[id] = false;
	} // end of while
	if(local_max_poly){
		pthread_mutex_lock(&output_lock);
		if(!max_poly){
			max_poly = local_max_poly->clone();
		}else if(local_max_poly->num_boundary_vertices()>max_poly->num_boundary_vertices()){
			delete max_poly;
			max_poly = local_max_poly;
		}
		pthread_mutex_unlock(&output_lock);
	}
	pthread_exit(NULL);
	return NULL;
}


void *process(void *args){
	int id = *(int *)args;
	pthread_mutex_lock(&output_lock);
	log("thread %d is started", id);
	pthread_mutex_unlock(&output_lock);

	vector<geometry_wrapper *> max_polies;
	PrecisionModel *pm = new PrecisionModel();
	GeometryFactory::unique_ptr gf = geos::geom::GeometryFactory::create(pm, OSM_SRID);
	WKTReader *wkt_reader = new WKTReader(gf.get());

	// Reading from the cache file
	/* Parsing polygon input */
	vector<string> result;

	while (!stop||is_working[id]) {
		if(!is_working[id]){
			usleep(10);
			continue;
		}
		tokenize(processing_line[id], result, "\t");
		try {
			Geometry *poly = wkt_reader->read(result[0]);
			geometry_wrapper *wrapper = new geometry_wrapper();
			wrapper->area = poly->getArea();
			wrapper->geom = poly;
			wrapper->raw = processing_line[id];
			update(max_polies,wrapper);
		} catch (const std::exception & ex ) {
			cerr << ex.what()<<endl;
		}
		result.clear();

		processing_line[id].clear();
		is_working[id] = false;
	} // end of while
	pthread_mutex_lock(&output_lock);
	for(geometry_wrapper *wrap:max_polies){
		update(global_top, wrap);
	}
	pthread_mutex_unlock(&output_lock);
	max_polies.clear();

	delete wkt_reader;
	delete pm;
	gf.release();
	pthread_exit(NULL);
	return NULL;
}


int main(int argc, char** argv) {
	if(argc>=3){
		k = atoi(argv[2]);
	}
	log("getting top %d",k);

	int num_threads = get_num_threads();
	if(argc>=2){
		num_threads = atoi(argv[1]);
	}
	pthread_t threads[num_threads];
	int id[num_threads];
	for(int i=0;i<num_threads;i++){
		id[i] = i;
		processing_line[i].clear();
		is_working[i] = false;
	}
	for(int i=0;i<num_threads;i++){
		pthread_create(&threads[i], NULL, process_wkt, (void *)&id[i]);
	}
	struct timeval start_time = get_cur_time();
	std::string input_line;
	int num_objects = 0;
	int next_report = 0;
	long processed_size = 0;

	while (getline(cin, input_line)) {
		while(true){
			bool assigned = false;
			for(int i=0;i<num_threads;i++){
				if(is_working[i]==false){
					processing_line[i] = input_line;
					is_working[i] = true;
					assigned = true;
					break;
				}
			}
			if(assigned){
				break;
			}
			usleep(10);
		}
		processed_size += input_line.size()+1;
		num_objects++;
		if(num_objects%1000000==0){
			log("processed %d objects", num_objects);
		}
	}
	stop = true;

	for(int i = 0; i < num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}

	logt("processed %d objects", start_time, num_objects);

	for(geometry_wrapper *geo:global_top){
		cout<<geo->geom->convexHull()->toString()<<endl;
		vector<string> strs;

		tokenize(geo->raw,strs,"\t");
		for(int i=1;i<strs.size();i++){
			if(i!=1){
				cout<<",";
			}
			cout<<strs[i];
		}
		cout<<endl<<endl;
		delete geo;
	}
	global_top.clear();

	if(max_poly){
		max_poly->print();
		cout<<max_poly->num_boundary_vertices()<<endl;
		delete max_poly;
	}

	pthread_exit(NULL);
	return true;
}
