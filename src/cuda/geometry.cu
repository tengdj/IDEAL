#include <cuda.h>
#include "../geometry/MyPolygon.h"
#include "mygpu.h"
#include "cuda_util.h"

// return the distance of two segments

__device__
bool check_contain(const double *polygon1, const double *polygon2, int num_vertices_1, int num_vertices_2){
	bool val = false;
	for(int p = 0;p<num_vertices_2-1;p++){
		double px = polygon2[p];
		double py = polygon2[num_vertices_2+p];
		for (int i = 0, j = 1; i < num_vertices_1-1; i++,j++) {
			// segment i->j intersect with line y=p.y
			double pyi = polygon1[num_vertices_1+i];
			double pyj = polygon1[num_vertices_1+j];
			if ((pyj>py) != (pyi>py))
			{
				double pxi = polygon1[i];
				double pxj = polygon1[j];
				double a = (pxj-pxi) / (pyj-pyi);
				if(px-pxi<a*(py-pyi)){
					val = !val;
				}
			}
		}
	}
	return val;
}

__global__
void contain_cuda(const double *poly1, const double *poly2, const uint *offset_size, int *ret, size_t num_pairs){

	// which polygon-polygon pair
	int pair_id = blockIdx.x*blockDim.x+threadIdx.x;
	if(pair_id>=num_pairs){
		return;
	}

	uint off1 = offset_size[pair_id*4];
	uint size1 = offset_size[pair_id*4+1];
	uint off2 = offset_size[pair_id*4+2];
	uint size2 = offset_size[pair_id*4+3];
//	printf("os: %d %d %d %d\n",off1, size1, off2,size2);
//
//	for(int i=0;i<size1;i++){
//		//printf("%f %f\n",(poly1+off1)[i],(poly1+off1)[i+size1]);
//	}
//	for(int i=0;i<size2;i++){
//		printf("%f %f\n",(poly2+off2)[i],(poly2+off2)[i+size1]);
//	}
	//printf("%d %d %d %d\n",off1, size1, off2, size2);
	ret[pair_id] = (int)check_contain(poly1+off1,poly2+off2,size1,size2);
}

/*
 * data: contains the segments of the meshes mentioned in this join.
 * offset_size:  contains the offset in the data for each batch, and the sizes of two data sets
 * result: for the returned results for each batch
 * batch_num: number of computed batches
 *
 * */
void contain_batch_gpu(gpu_info *gpu, double *data, uint *offset_size, int *result, size_t total_vertice_num, int pair_num){

	assert(gpu);
	cudaSetDevice(gpu->device_id);
	struct timeval start = get_cur_time();

	// space for the results in GPU
	int *d_ret = gpu->get_result(sizeof(int)*pair_num);
	// space for the offset and size information in GPU
	uint *d_os = gpu->get_os(sizeof(uint)*pair_num*4);
	double *d_poly1 = gpu->source_data;
	double *d_poly2 = gpu->get_data(total_vertice_num*2*sizeof(double));

	CUDA_SAFE_CALL(cudaMemcpy(d_poly2, data, total_vertice_num*2*sizeof(double), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_os, offset_size, pair_num*4*sizeof(uint), cudaMemcpyHostToDevice));
	//logt("allocating data", start);

	// compute the vectors of segments in data, save to d_vec
	contain_cuda<<<pair_num/1024+1,1024>>>(d_poly1, d_poly2, d_os, d_ret, pair_num);
	check_execution();
	cudaDeviceSynchronize();
	//logt("computations", start);
	CUDA_SAFE_CALL(cudaMemcpy(result, d_ret, pair_num*sizeof(int), cudaMemcpyDeviceToHost));
	//logt("copy data out", start);
}

void load_source_togpu(gpu_info *gpu, vector<MyPolygon *> &source){
	size_t source_size = 0;
	for(MyPolygon *p:source){
		source_size += 2*p->boundary->num_vertices*sizeof(double);
	}
	gpu->get_source(source_size);
	source_size = 0;
	for(MyPolygon *p:source){
		int num_vertices = p->boundary->num_vertices;
		p->offset = source_size/sizeof(double);
		CUDA_SAFE_CALL(cudaMemcpy((char *)(gpu->source_data)+source_size, (char *)(p->boundary->x), num_vertices*sizeof(double),cudaMemcpyHostToDevice));
		source_size += num_vertices*sizeof(double);
		CUDA_SAFE_CALL(cudaMemcpy((char *)(gpu->source_data)+source_size, (char *)(p->boundary->y), num_vertices*sizeof(double),cudaMemcpyHostToDevice));
		source_size += num_vertices*sizeof(double);
	}
}


