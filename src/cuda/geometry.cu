
#include "cuda_util.cuh"
#include "MyPolygon.h"

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
	ret[pair_id] = (int)check_contain(poly1+off1,poly2+off2,size1,size2);
}

///*
// * data: contains the segments of the meshes mentioned in this join.
// * offset_size:  contains the offset in the data for each batch, and the sizes of two data sets
// * result: for the returned results for each batch
// * batch_num: number of computed batches
// *
// * */
//void contain_batch_gpu(gpu_info *gpu, double *data, uint *offset_size, int *result, size_t total_vertice_num, int pair_num){
//
//	assert(gpu);
//	cudaSetDevice(gpu->device_id);
//	struct timeval start = get_cur_time();
//
//	// space for the results in GPU
//	int *d_ret = gpu->get_result(sizeof(int)*pair_num);
//	// space for the offset and size information in GPU
//	uint *d_os = gpu->get_os(sizeof(uint)*pair_num*4);
//	double *d_poly1 = gpu->source_data;
//	double *d_poly2 = gpu->get_data(total_vertice_num*2*sizeof(double));
//
//	CUDA_SAFE_CALL(cudaMemcpy(d_poly2, data, total_vertice_num*2*sizeof(double), cudaMemcpyHostToDevice));
//	CUDA_SAFE_CALL(cudaMemcpy(d_os, offset_size, pair_num*4*sizeof(uint), cudaMemcpyHostToDevice));
//	//logt("allocating data", start);
//
//	// compute the vectors of segments in data, save to d_vec
//	contain_cuda<<<pair_num/1024+1,1024>>>(d_poly1, d_poly2, d_os, d_ret, pair_num);
//	check_execution();
//	cudaDeviceSynchronize();
//	//logt("computations", start);
//	CUDA_SAFE_CALL(cudaMemcpy(result, d_ret, pair_num*sizeof(int), cudaMemcpyDeviceToHost));
//	//logt("copy data out", start);
//}

__device__
double cuda_degree_per_kilometer_longitude(double latitude, double *degree_per_kilometer){
	double absla = abs(latitude);
	assert(absla<=90);
	if(absla==90){
		absla = 89.9;
	}
	return degree_per_kilometer[(int)(absla*10)];
}


__device__
double cuda_point_to_segment_distance(const Point &p, const Point &p1, const Point &p2, double *degree_per_kilometer) {

  double A = p.x - p1.x;
  double B = p.y - p1.y;
  double C = p2.x - p1.x;
  double D = p2.y - p1.y;

  double dot = A * C + B * D;
  double len_sq = C * C + D * D;
  double param = -1;
  if (len_sq != 0) //in case of 0 length line
      param = dot / len_sq;

  double xx, yy;

  if (param < 0) {
    xx = p1.x;
    yy = p1.y;
  } else if (param > 1) {
    xx = p2.x;
    yy = p2.y;
  } else {
    xx = p1.x + param * C;
    yy = p1.y + param * D;
  }

  double dx = p.x - xx;
  double dy = p.y - yy;
  dx = dx/cuda_degree_per_kilometer_longitude(p.y, degree_per_kilometer);
  dy = dy/degree_per_kilometer_latitude;

  return sqrt(dx * dx + dy * dy);
}

__device__
void atomicMin_double(double* address, double val)
{
    unsigned long long int* address_as_ull = (unsigned long long int*) address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
            __double_as_longlong(fmin(val, __longlong_as_double(assumed))));
    } while (assumed != old);
}

__global__
void cuda_distance(double *dist, const Point p, const Point *vs, size_t vs_length,double *degree_per_kilometer){

	// which polygon-polygon pair
	int pair_id = blockIdx.x*blockDim.x+threadIdx.x;
	if(pair_id>=vs_length){
		return;
	}
	double d = cuda_point_to_segment_distance(p, vs[pair_id], vs[pair_id+1],degree_per_kilometer);
	//if(d>0.00001)
	{
		atomicMin_double(dist, d);
	}
}

double point_to_segment_sequence_distance_gpu(Point &p, Point *vs, size_t seq_len, bool geography){
	assert(gpus.size()>0);
	gpu_info *gpu = gpus[0];
	Point *vs_d = (Point *)gpu->allocate(seq_len* sizeof(Point));
	double *dist_d = (double *)gpu->allocate(sizeof(double));
	double dist = DBL_MAX;
	CUDA_SAFE_CALL(cudaMemcpy((char *)dist_d, (char *)&dist, sizeof(double),cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy((char *)vs_d, (char *)vs, seq_len*sizeof(Point),cudaMemcpyHostToDevice));
	cuda_distance<<<seq_len/1024+1, 1024>>>(dist_d, p, vs_d, seq_len,gpu->degree_per_kilometer);
	CUDA_SAFE_CALL(cudaMemcpy(&dist, dist_d, sizeof(double), cudaMemcpyDeviceToHost));
	gpu->free((char *)dist_d, sizeof(double));
	gpu->free((char *)vs_d, seq_len* sizeof(Point));
	return dist;
}



//void load_source_togpu(gpu_info *gpu, vector<MyPolygon *> &source){
//	size_t source_size = 0;
//	for(MyPolygon *p:source){
//		source_size += 2*p->boundary->num_vertices*sizeof(double);
//	}
//	gpu->get_source(source_size);
//	source_size = 0;
//	for(MyPolygon *p:source){
//		int num_vertices = p->boundary->num_vertices;
//		p->offset = source_size/sizeof(double);
//		CUDA_SAFE_CALL(cudaMemcpy((char *)(gpu->source_data)+source_size, (char *)(p->boundary->x), num_vertices*sizeof(double),cudaMemcpyHostToDevice));
//		source_size += num_vertices*sizeof(double);
//		CUDA_SAFE_CALL(cudaMemcpy((char *)(gpu->source_data)+source_size, (char *)(p->boundary->y), num_vertices*sizeof(double),cudaMemcpyHostToDevice));
//		source_size += num_vertices*sizeof(double);
//	}
//}
////
////
////void contain_batch_gpu(gpu_info *gpu, double *data, uint *offset_size, int *result, size_t total_vertice_num, int pair_num);
//////void load_source_togpu(gpu_info *gpu, vector<MyPolygon *> &source);
////
//int process_with_gpu(gpu_info *gpu, query_context *ctx){
//
//	int pair_num = ctx->candidates.size();
//	if(pair_num==0){
//		return 0;
//	}
//	int *result = new int[pair_num];
//	uint *offset_size = new uint[pair_num*4];
//	uint dataoffset = 0;
//	uint total_num_vertices = 0;
//	for(int i=0;i<pair_num;i++){
//		if(i==0||ctx->candidates[i].second->getid()!=ctx->candidates[i-1].second->getid()){
//			total_num_vertices += ctx->candidates[i].second->boundary->num_vertices;
//		}
//	}
//
//	double *tmpdata = new double[total_num_vertices*2];
//	for(int i=0;i<pair_num;i++){
//		offset_size[i*4] = ctx->candidates[i].first->offset;
//		offset_size[i*4+1] = ctx->candidates[i].first->boundary->num_vertices;
//		offset_size[i*4+3] = ctx->candidates[i].second->boundary->num_vertices;
//		if(i==0||ctx->candidates[i].second->getid()!=ctx->candidates[i-1].second->getid()){
//			offset_size[i*4+2] = dataoffset;
//			int num_vertices = ctx->candidates[i].second->boundary->num_vertices;
//			memcpy((char *)(tmpdata+dataoffset), (char *)(ctx->candidates[i].second->boundary->x), num_vertices*sizeof(double));
//			dataoffset += num_vertices;
//			memcpy((char *)(tmpdata+dataoffset), (char *)(ctx->candidates[i].second->boundary->y), num_vertices*sizeof(double));
//			dataoffset += num_vertices;
//		}else{
//			offset_size[i*4+2] = dataoffset-offset_size[i*4+3]*2;
//		}
//	}
//	assert(dataoffset==total_num_vertices*2);
//	ctx->candidates.clear();
//
//	int found = 0;
//	pthread_mutex_lock(&gpu->lock);
//	contain_batch_gpu(gpu,tmpdata,offset_size,result,total_num_vertices,pair_num);
//	pthread_mutex_unlock(&gpu->lock);
//	for(int i=0;i<pair_num;i++){
//		found += result[i];
//	}
//	delete []result;
//	delete []offset_size;
//	delete []tmpdata;
//	return found;
//}

