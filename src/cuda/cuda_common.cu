/*
 *
 * with some common gpu related operations
 * */

#include <cuda.h>
#include "mygpu.h"
#include "cuda_util.h"
#include "../include/util.h"

using namespace std;



vector<gpu_info *> get_gpus(){
	vector<gpu_info *> gpus;
	int num_gpus = 0;
	cudaGetDeviceCount(&num_gpus);
	for (int i = 0; i < num_gpus; i++) {
		cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop, i);
		gpu_info *info = new gpu_info();
		info->busy = false;
		info->mem_size = prop.totalGlobalMem/1024/1024*4/5;
		info->device_id = i;
		// we allocate 2G mem for each gpu
//		if(info->mem_size>2048){
//			info->mem_size = 2048;
//		}
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


void gpu_info::init(){
}

double *gpu_info::get_source(size_t ss){
	cudaSetDevice(this->device_id);
	if(this->source_data&&this->source_size<ss){
		CUDA_SAFE_CALL(cudaFree(this->source_data));
		this->source_size = 0;
		this->source_data = NULL;
	}
	if(!this->source_data){
		CUDA_SAFE_CALL(cudaMalloc((void **)&source_data, ss));
		assert(this->source_data);
		this->source_size = ss;
	}
	return this->source_data;
}
double *gpu_info::get_data(size_t ds){
	cudaSetDevice(this->device_id);
	if(this->d_data&&this->data_size<ds){
		CUDA_SAFE_CALL(cudaFree(this->d_data));
		this->data_size = 0;
		this->d_data = NULL;
	}
	if(!this->d_data){
		CUDA_SAFE_CALL(cudaMalloc((void **)&this->d_data, ds));
		assert(this->d_data);
		this->data_size = ds;
	}
	return this->d_data;
}
int *gpu_info::get_result(size_t rs){
	cudaSetDevice(this->device_id);
	if(this->result&&this->result_size<rs){
		CUDA_SAFE_CALL(cudaFree(this->result));
		this->result_size = 0;
		this->result = NULL;
	}
	if(!this->result){
		CUDA_SAFE_CALL(cudaMalloc((void **)&this->result, rs));
		assert(this->result);
		this->result_size = rs;
	}
	return this->result;
}
uint *gpu_info::get_os(size_t os){
	cudaSetDevice(this->device_id);
	if(this->offset_size&&this->os_size<os){
		CUDA_SAFE_CALL(cudaFree(this->offset_size));
		this->os_size = 0;
		this->offset_size = NULL;
	}
	if(!this->offset_size){
		CUDA_SAFE_CALL(cudaMalloc((void **)&this->offset_size, os));
		assert(this->offset_size);
		this->os_size = os;
	}
	return this->offset_size;
}



gpu_info::~gpu_info(){
	cudaSetDevice(this->device_id);
	if(this->d_data){
		CUDA_SAFE_CALL(cudaFree(this->d_data));
		this->d_data = NULL;
	}
	if(this->source_data){
		CUDA_SAFE_CALL(cudaFree(this->source_data));
		this->source_data = NULL;
	}
	if(this->result){
		CUDA_SAFE_CALL(cudaFree(this->result));
		this->result = NULL;
	}
	if(this->offset_size){
		CUDA_SAFE_CALL(cudaFree(this->offset_size));
		this->offset_size = NULL;
	}
}

