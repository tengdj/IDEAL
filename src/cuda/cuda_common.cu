/*
 *
 * with some common gpu related operations
 * */

#include <cuda.h>
#include "cuda_util.cuh"
#include "util.h"

using namespace std;

vector<gpu_info *> gpus;

void init_gpus(){
	assert(gpus.size()==0&&"gpu need to be initialized once");
	gpus = get_gpus();
}

void release_gpus(){
	for(gpu_info *g:gpus){
		delete g;
	}
	gpus.clear();
}

vector<gpu_info *> get_gpus(){
	vector<gpu_info *> gpus;
	int num_gpus = 0;
	cudaGetDeviceCount(&num_gpus);
	for (int i = 0; i < num_gpus; i++) {
		cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop, i);
		gpu_info *info = new gpu_info(i, (size_t )prop.totalGlobalMem*4/5);
		gpus.push_back(info);
	}
	return gpus;
}

void print_gpus(){
	int num_gpus = 0;
	cudaGetDeviceCount(&num_gpus);
	for (int i = 0; i < num_gpus; i++) {
		cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop, i);
		log("Device Number: %d", i);
		log("  Device name: %s", prop.name);
		log("  Memory Clock Rate (KHz): %d", prop.memoryClockRate);
		log("  Memory Bus Width (bits): %d", prop.memoryBusWidth);
		log("  Peak Memory Bandwidth (GB/s): %f",
				2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
		log("  Memory size (MB): %ld\n", prop.totalGlobalMem/1024/1024);
	}
}

char *gpu_info::allocate(size_t ss){

	lock();
	cudaSetDevice(device_id);
	size_t offset = 0;
	// in odd loop
	if(allocated_end>=allocated_begin){
		if(allocated_end+ss<mem_size){
			offset = allocated_end;
			allocated_end += ss;
		}else if(allocated_begin>=ss){
			// the space between allocated_end to the end is wasted.
			freed[allocated_end] = 0;
			offset = 0;
			allocated_end = ss;
		}else{
			assert(false && "out of memory");
		}
	}else{
		// in even loop, end in front of begin
		if(allocated_end+ss<allocated_begin){
			offset = allocated_end;
			allocated_end += ss;
		}else{
			assert(false && "out of memory");
		}
	}
	unlock();
	return d_data+offset;
}

void gpu_info::free(char *daddr, size_t sz){
	size_t free_begin = (size_t)(daddr-d_data);
	size_t free_end = free_begin+sz;
	//log("freeing %ld to %ld", free_begin, free_end);
	lock();
	cudaSetDevice(this->device_id);
	if(free_begin==allocated_begin){
		do{
			allocated_begin = free_end;
			// some neighbor memory space on the right size already freed, move forward
			if(freed.find(free_end)!=freed.end()){
				size_t new_end = freed[free_end];
				freed.erase(free_end);
				free_end = new_end;
			}else{
				break;
			}
		}while(true);
	}else{
		// some previous memory space is not freed yet
		freed[free_begin] = free_end;
	}
	if(allocated_begin==allocated_end&&allocated_begin!=0){
		allocated_begin = 0;
		allocated_end = 0;
	}
	unlock();
	//log("after freeing %ld to %ld", allocated_begin, allocated_end);
}

gpu_info::gpu_info(int gid, size_t memsize){
	device_id = gid;
	mem_size = memsize;
	busy = false;
	pthread_mutex_init(&lk,NULL);
	CUDA_SAFE_CALL(cudaMalloc(&d_data, mem_size));
	size_t dgr_size = sizeof(degree_per_kilometer_longitude_arr);
	CUDA_SAFE_CALL(cudaMalloc(&degree_per_kilometer, dgr_size));
	CUDA_SAFE_CALL(cudaMemcpy((char *)degree_per_kilometer, (char *)degree_per_kilometer_longitude_arr,dgr_size,cudaMemcpyHostToDevice));
}

gpu_info::~gpu_info(){
	cudaSetDevice(this->device_id);
	if(this->d_data){
		CUDA_SAFE_CALL(cudaFree(this->d_data));
		this->d_data = NULL;
	}
	if(this->degree_per_kilometer){
		CUDA_SAFE_CALL(cudaFree(this->degree_per_kilometer));
		this->degree_per_kilometer = NULL;
	}
}

void gpu_info::lock(){
	pthread_mutex_lock(&lk);
}

void gpu_info::unlock(){
	pthread_mutex_unlock(&lk);
}

inline double cuda_degree_per_kilometer_longitude_calculate(double latitude){
	double absla = abs(latitude);
	assert(absla<=90);
	if(absla==90){
		absla = 89;
	}
	return 360.0/(sin((90-absla)*PI/180)*40076);
}

inline double cuda_degree_per_kilometer_longitude(double latitude){
	double absla = abs(latitude);
	assert(absla<=90);
	if(absla==90){
		absla = 89.9;
	}
	return degree_per_kilometer_longitude_arr[(int)(absla*10)];
}

