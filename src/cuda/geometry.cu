#include "cuda_util.h"
#include "Ideal.h"
#include "mygpu.h"

#define BLOCK_SIZE 1024

struct Threadwork{
	int source_pixid = 0;
	int target_pixid = 0;
	int pair_id = 0;
};

__device__ int gpu_get_id(int x, int y, int dimx){
	return y * (dimx+1) + x;
}

__device__ int  gpu_double_to_int(double val){
	int vi = (int)val;
	if(abs(1.0*(vi+1)-val)<0.00000001){
		vi++;
	}
	return vi;
}

__device__ int get_offset_x(double s_xval, double t_xval, double step_x, int dimx){
	int x = gpu_double_to_int((t_xval-s_xval)/step_x);
	return min(max(x, 0), dimx);
}

__device__ int get_offset_y(double s_yval, double t_yval, double step_y, int dimy){
	int y = gpu_double_to_int((t_yval-s_yval)/step_y);
	return min(max(y, 0), dimy);
}

__device__ PartitionStatus gpu_show_status(uint8_t *status, uint start, int id){
	uint8_t st = (status+start)[id / 4];
	int pos = id % 4 * 2;   // The multiplication by 2 is because each status occupies 2 bits.	
	st &= ((uint8_t)3 << pos);
	st >>= pos;
	if(st == 0) return OUT;
	if(st == 3) return IN;
	return BORDER;
}

// __device__ box gpu_get_pixel_box(box &mbr, int x, int y, double step_x, double step_y){
//     box bx;
//     const double start_x = mbr.low[0];
//     const double start_y = mbr.low[1];

//     double lowx = start_x + x * step_x;
//     double lowy = start_y + y * step_y;
//     double highx = start_x + (x + 1) * step_x;
//     double highy = start_y + (y + 1) * step_y;

//     bx.low[0] = lowx;
//     bx.low[1] = lowy;
//     bx.high[0] = highx;
//     bx.high[1] = highy;
//     return bx;
// }

// __global__ void kernel_retrieve_pixels(pair<IdealOffset, IdealOffset> *d_pairs, Idealinfo *d_info, uint8_t *d_status, uint size, uint8_t *d_pixels, uint *d_pixels_idx, uint *d_pixels_offset){
// 	const int x = blockIdx.x * blockDim.x + threadIdx.x;
// 	if(x < size){
// 		pair<IdealOffset, IdealOffset> temp_pair = d_pairs[x];
// 		IdealOffset source = temp_pair.first;
// 		IdealOffset target = temp_pair.second;	
// 		box s_mbr = d_info[source.info_start].mbr, t_mbr = d_info[target.info_start].mbr;
// 		double step_x = d_info[source.info_start].step_x, step_y = d_info[source.info_start].step_y;
// 		int dimx = d_info[source.info_start].dimx, dimy = d_info[source.info_start].dimy;
// 		// printf("source MBR: %lf %lf %lf %lf\n", s_mbr.low[0], s_mbr.low[1], s_mbr.high[0], s_mbr.high[1]);
// 		// printf("target MBR: %lf %lf %lf %lf\n", t_mbr.low[0], t_mbr.low[1], t_mbr.high[0], t_mbr.high[1]);
// 		// printf("%d %d\n", dimx, dimy);
// 		// printf("%lf, %lf\n", step_x, step_y); 
// 		int start_x = get_offset_x(s_mbr.low[0], t_mbr.low[0], step_x, dimx);
// 		int start_y = get_offset_y(s_mbr.low[1], t_mbr.low[1], step_y, dimy);
// 		int end_x = get_offset_x(s_mbr.high[0], t_mbr.high[0], step_x, dimx);
// 		int end_y = get_offset_y(s_mbr.high[1], t_mbr.high[1], step_y, dimy);
		
// 		uint idx = atomicAdd(d_pixels_idx, (end_x-start_x+1)*(end_y-start_y+1));
// 		for(int i=start_x;i<=end_x;i++){
// 			for(int j=start_y;j<=end_y;j++){
// 				int id = gpu_get_id(i , j, dimx);
// 				d_pixels[idx ++] = gpu_show_status(d_status, source.status_start, id);
// 				d_pixels_offset[x] = idx;
// 			}
// 		}
// 		// printf("%d FROM GPU\n", *d_pixels_idx);
// 		// atomicAdd(d_pixels_idx, 1);
// 		// printf("%d FROM GPU\n", x);
// 	}

// }

__global__ void kernel_filter(pair<IdealOffset, IdealOffset> *d_pairs,Idealinfo *d_info, uint8_t *d_status, uint size, uint8_t *resultmap){
	const int x = blockIdx.x * blockDim.x + threadIdx.x;
	if(x < size){
		pair<IdealOffset, IdealOffset> temp_pair = d_pairs[x];
		IdealOffset source = temp_pair.first;
		IdealOffset target = temp_pair.second;		
		
		box s_mbr = d_info[source.info_start].mbr, t_mbr = d_info[target.info_start].mbr;
		double s_step_x = d_info[source.info_start].step_x, s_step_y = d_info[source.info_start].step_y;
		int s_dimx = d_info[source.info_start].dimx, s_dimy = d_info[source.info_start].dimy;

		int start_x = get_offset_x(s_mbr.low[0], t_mbr.low[0], s_step_x, s_dimx);
		int start_y = get_offset_y(s_mbr.low[1], t_mbr.low[1], s_step_y, s_dimy);
		int end_x = get_offset_x(s_mbr.low[0], t_mbr.high[0], s_step_x, s_dimx);
		int end_y = get_offset_y(s_mbr.low[1], t_mbr.high[1], s_step_y, s_dimy);

		uint itn = 0, etn = 0;
		for(int i=start_x;i<=end_x;i++){
			for(int j=start_y;j<=end_y;j++){
				int id = gpu_get_id(i , j, s_dimx);
				if(gpu_show_status(d_status, source.status_start, id) == IN) itn ++;
				else if(gpu_show_status(d_status, source.status_start, id) == OUT) etn ++;
			}
		}
		if(itn == (end_x-start_x+1)*(end_y-start_y+1)){
			resultmap[x] = 1;
		}
		if(etn == (end_x-start_x+1)*(end_y-start_y+1)){
			resultmap[x] = 2;
		}
	}
}

// __device__ int gpu_sgn(const double& x) {
// 	return x >= 0 ? x ? 1 : 0 : -1;
// }

__device__ int gpu_sgn(const double& x) {
    return (x > 0) - (x < 0);
}




__device__ bool gpu_inter1(double a, double b, double c, double d) {

	double tmp;
    if (a > b){
    	tmp = a;
    	a = b;
    	b = tmp;
    }
    if (c > d){
    	tmp = c;
    	c = d;
    	d = tmp;
    }
    return max(a, c) <= min(b, d);
}

// __device__ double gpu_cross(const Point& p, const Point& q) {
// 	return q.x * p.y - q.y * p.x;
// }


__device__ double gpu_cross(const Point &a, const Point &b, const Point& c){
	double x1, y1, x2, y2;
	x1 = a.x - c.x;
	y1 = a.y - c.y;
	x2 = b.x - c.x;
	y2 = b.y - c.y;
	return x1 * y2 - y1 * x2;
}

__device__ bool gpu_segment_intersect(const Point& a, const Point& b, const Point& c, const Point& d) {
    if (gpu_cross(a, d, c) == 0 && gpu_cross(b, d, c) == 0)
        return gpu_inter1(a.x, b.x, c.x, d.x) && gpu_inter1(a.y, b.y, c.y, d.y);

    return gpu_sgn(gpu_cross(b, c, a)) != gpu_sgn(gpu_cross(b, d, a)) &&
          gpu_sgn(gpu_cross(d, a, c)) != gpu_sgn(gpu_cross(d, b, c));

}



__device__ bool gpu_segment_intersect_batch(Point *p1, Point *p2, int s1, int s2){
	for(int i=0;i<s1;i++){
		for(int j=0;j<s2;j++){
			if(gpu_segment_intersect(p1[i],p1[(i+1)%s1],p2[j],p2[(j+1)%s2])){
				return true;
			}
		}
	}
	return false;
}

__global__ void kernel_refinement(Threadwork *d_threadwork, pair<IdealOffset, IdealOffset> *d_pairs, uint16_t *d_offset, EdgeSeq *d_edge_sequences, Point *d_vertices, uint size, uint8_t* resultmap){
	const int x = blockIdx.x * blockDim.x + threadIdx.x;
	if(x < size){
		int p = d_threadwork[x].source_pixid;
		int p2 = d_threadwork[x].target_pixid;
		int pair_id = d_threadwork[x].pair_id;
		if(resultmap[pair_id] != 0) return;

		pair<IdealOffset, IdealOffset> temp_pair = d_pairs[pair_id];
		IdealOffset source = temp_pair.first;
		IdealOffset target = temp_pair.second;	

		uint s_offset_start = source.offset_start, t_offset_start = target.offset_start; 
		uint s_edge_sequences_start = source.edge_sequences_start, t_edge_sequences_start = target.edge_sequences_start; 
		int s_num_sequence = (d_offset+s_offset_start)[p + 1] - (d_offset+s_offset_start)[p];
		int t_num_sequence = (d_offset+t_offset_start)[p2 + 1] - (d_offset+t_offset_start)[p2];
		uint s_vertices_start = source.vertices_start, t_vertices_start = target.vertices_start;
		// int idx = atomicAdd(d_segment_size, 1U*s_num_sequence*t_num_sequence);
		for(int i = 0; i < s_num_sequence; ++ i){
			EdgeSeq r = (d_edge_sequences+s_edge_sequences_start)[(d_offset+s_offset_start)[p] + i]; 
			for(int j = 0; j < t_num_sequence; ++ j){
				EdgeSeq r2 = (d_edge_sequences+t_edge_sequences_start)[(d_offset+t_offset_start)[p2] + j];
				if(gpu_segment_intersect_batch((d_vertices+s_vertices_start+r.start), (d_vertices+t_vertices_start+r2.start), r.length, r2.length)){
					resultmap[pair_id] = 3;
					return;
				}
			}
		}
	}
}

void initialize_pairs(query_context *gctx, pair<IdealOffset, IdealOffset> *h_pairs, uint size, std::atomic<int> &atomic_i) {
    while (true) {
        int i = atomic_i.fetch_add(1);
        if (i >= size) break;
        
        Ideal *source = gctx->ideal_pairs[i].first;
        Ideal *target = gctx->ideal_pairs[i].second;
        h_pairs[i] = {*source->idealoffset, *target->idealoffset};
    }
}

void process_pixels(query_context *gctx, uint8_t *h_resultmap, Threadwork *h_threadwork, std::atomic<int> &atomic_id, int start, int end) {
    for (int i = start; i < end; ++i) {
        if (h_resultmap[i] == 0) {
            Ideal *source = gctx->ideal_pairs[i].first;
            Ideal *target = gctx->ideal_pairs[i].second;
            vector<int> pxs = source->retrieve_pixels(target->getMBB());
            for (auto p : pxs) {
                box bx = source->get_pixel_box(source->get_x(p), source->get_y(p));
                vector<int> tpxs = target->retrieve_pixels(&bx);
                for (auto p2 : tpxs) {
                    if (source->show_status(p) == BORDER && target->show_status(p2) == BORDER) {
                        int id = atomic_id.fetch_add(1);
                        h_threadwork[id] = {p, p2, i};
                    }
                }
            }
        }
    }
}

uint cuda_contain(query_context *gctx){
	uint size = gctx->ideal_pairs.size();
	
	pair<IdealOffset, IdealOffset> *h_pairs = nullptr;
	pair<IdealOffset, IdealOffset> *d_pairs = nullptr;

	h_pairs = new pair<IdealOffset, IdealOffset>[size];
	
	for(int i = 0; i < size; ++ i){
		Ideal *source = gctx->ideal_pairs[i].first;
		Ideal *target = gctx->ideal_pairs[i].second;
		h_pairs[i] = {*source->idealoffset, *target->idealoffset};
	}

	CUDA_SAFE_CALL(cudaMalloc((void**) &d_pairs, size * sizeof(pair<IdealOffset, IdealOffset>)));
	CUDA_SAFE_CALL(cudaMemcpy(d_pairs, h_pairs, size *  sizeof(pair<IdealOffset, IdealOffset>), cudaMemcpyHostToDevice));

	uint8_t *d_resultmap = nullptr;
	CUDA_SAFE_CALL(cudaMalloc((void **) &d_resultmap, size * sizeof(uint8_t)));
	CUDA_SAFE_CALL(cudaMemset(d_resultmap, 0, size * sizeof(uint8_t)));

	int grid_size_x = (size + BLOCK_SIZE - 1) / BLOCK_SIZE;
	dim3 block_size(BLOCK_SIZE, 1, 1);
	dim3 grid_size(grid_size_x, 1, 1);
	kernel_filter<<<grid_size, block_size>>>(d_pairs, gctx->d_info, gctx->d_status, size, d_resultmap);
	cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess) {
        printf("FILTER CUDA error: %s\n", cudaGetErrorString(error));
    }

	cudaDeviceSynchronize();

	uint8_t *h_resultmap = new uint8_t[size];
	CUDA_SAFE_CALL(cudaMemcpy(h_resultmap, d_resultmap, size * sizeof(uint8_t), cudaMemcpyDeviceToHost));

	Threadwork *h_threadwork = new Threadwork[8*1024*1024];
	Threadwork *d_threadwork = nullptr;
	CUDA_SAFE_CALL(cudaMalloc((void **) &d_threadwork, 8*1024*1024*sizeof(Threadwork)));

	std::atomic<int> atomic_id(0);

	int num_threads = std::thread::hardware_concurrency();
    std::vector<std::thread> threads;
    int chunk_size = (size + num_threads - 1) / num_threads;
    for (int i = 0; i < num_threads; ++i) {
        int start = i * chunk_size;
        int end = std::min((i + 1) * chunk_size, static_cast<int>(size));
        threads.emplace_back(process_pixels, gctx, h_resultmap, h_threadwork, std::ref(atomic_id), start, end);
    }

    for (auto &thread : threads) {
        thread.join();
    }

	int id = atomic_id.load();

	CUDA_SAFE_CALL(cudaMemcpy(d_threadwork, h_threadwork, 8*1024*1024*sizeof(Threadwork), cudaMemcpyHostToDevice));

	grid_size_x = (id + BLOCK_SIZE - 1) / BLOCK_SIZE;
	grid_size.x = grid_size_x;
	
	kernel_refinement<<<grid_size, block_size>>>(d_threadwork, d_pairs, gctx->d_offset, gctx->d_edge_sequences, gctx->d_vertices, id, d_resultmap);

	error = cudaGetLastError();
    if (error != cudaSuccess) {
        printf("REFINEMENT CUDA error: %s\n", cudaGetErrorString(error));
    }

	cudaDeviceSynchronize();

	CUDA_SAFE_CALL(cudaMemcpy(h_resultmap, d_resultmap, size * sizeof(uint8_t), cudaMemcpyDeviceToHost));

	int found = 0;
	for(int i = 0; i < size; ++ i ){
		if(h_resultmap[i] == 1 || h_resultmap[i] == 0) found ++;
	}

	delete []h_pairs;
	delete []h_resultmap;
	delete []h_threadwork;

	CUDA_SAFE_CALL(cudaFree(d_pairs));
	CUDA_SAFE_CALL(cudaFree(d_resultmap));
	CUDA_SAFE_CALL(cudaFree(d_threadwork));

	return found;
}

// __device__
// bool check_contain(const double *polygon1, const double *polygon2, int num_vertices_1, int num_vertices_2){
// 	bool val = false;
// 	for(int p = 0;p<num_vertices_2-1;p++){
// 		double px = polygon2[p];
// 		double py = polygon2[num_vertices_2+p];
// 		for (int i = 0, j = 1; i < num_vertices_1-1; i++,j++) {
// 			// segment i->j intersect with line y=p.y
// 			double pyi = polygon1[num_vertices_1+i];
// 			double pyj = polygon1[num_vertices_1+j];
// 			if ((pyj>py) != (pyi>py))
// 			{
// 				double pxi = polygon1[i];
// 				double pxj = polygon1[j];
// 				double a = (pxj-pxi) / (pyj-pyi);
// 				if(px-pxi<a*(py-pyi)){
// 					val = !val;
// 				}
// 			}
// 		}
// 	}
// 	return val;
// }

// __global__
// void contain_cuda(const double *poly1, const double *poly2, const uint *offset_size, int *ret, size_t num_pairs){

// 	// which polygon-polygon pair
// 	int pair_id = blockIdx.x*blockDim.x+threadIdx.x;
// 	if(pair_id>=num_pairs){
// 		return;
// 	}

// 	uint off1 = offset_size[pair_id*4];
// 	uint size1 = offset_size[pair_id*4+1];
// 	uint off2 = offset_size[pair_id*4+2];
// 	uint size2 = offset_size[pair_id*4+3];
// 	ret[pair_id] = (int)check_contain(poly1+off1,poly2+off2,size1,size2);
// }

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

// __device__
// double cuda_degree_per_kilometer_longitude(double latitude, double *degree_per_kilometer){
// 	double absla = abs(latitude);
// 	assert(absla<=90);
// 	if(absla==90){
// 		absla = 89.9;
// 	}
// 	return degree_per_kilometer[(int)(absla*10)];
// }


// __device__
// double cuda_point_to_segment_distance(const Point &p, const Point &p1, const Point &p2, double *degree_per_kilometer) {

//   double A = p.x - p1.x;
//   double B = p.y - p1.y;
//   double C = p2.x - p1.x;
//   double D = p2.y - p1.y;

//   double dot = A * C + B * D;
//   double len_sq = C * C + D * D;
//   double param = -1;
//   if (len_sq != 0) //in case of 0 length line
//       param = dot / len_sq;

//   double xx, yy;

//   if (param < 0) {
//     xx = p1.x;
//     yy = p1.y;
//   } else if (param > 1) {
//     xx = p2.x;
//     yy = p2.y;
//   } else {
//     xx = p1.x + param * C;
//     yy = p1.y + param * D;
//   }

//   double dx = p.x - xx;
//   double dy = p.y - yy;
//   dx = dx/cuda_degree_per_kilometer_longitude(p.y, degree_per_kilometer);
//   dy = dy/degree_per_kilometer_latitude;

//   return sqrt(dx * dx + dy * dy);
// }

// __device__
// void atomicMin_double(double* address, double val)
// {
//     unsigned long long int* address_as_ull = (unsigned long long int*) address;
//     unsigned long long int old = *address_as_ull, assumed;
//     do {
//         assumed = old;
//         old = atomicCAS(address_as_ull, assumed,
//             __double_as_longlong(fmin(val, __longlong_as_double(assumed))));
//     } while (assumed != old);
// }

// __global__
// void cuda_distance(double *dist, const Point p, const Point *vs, size_t vs_length,double *degree_per_kilometer){

// 	// which polygon-polygon pair
// 	int pair_id = blockIdx.x*blockDim.x+threadIdx.x;
// 	if(pair_id>=vs_length){
// 		return;
// 	}
// 	double d = cuda_point_to_segment_distance(p, vs[pair_id], vs[pair_id+1],degree_per_kilometer);
// 	//if(d>0.00001)
// 	{
// 		atomicMin_double(dist, d);
// 	}
// }

// double point_to_segment_sequence_distance_gpu(Point &p, Point *vs, size_t seq_len, bool geography){
// 	assert(gpus.size()>0);
// 	gpu_info *gpu = gpus[0];
// 	Point *vs_d = (Point *)gpu->allocate(seq_len* sizeof(Point));
// 	double *dist_d = (double *)gpu->allocate(sizeof(double));
// 	double dist = DBL_MAX;
// 	CUDA_SAFE_CALL(cudaMemcpy((char *)dist_d, (char *)&dist, sizeof(double),cudaMemcpyHostToDevice));
// 	CUDA_SAFE_CALL(cudaMemcpy((char *)vs_d, (char *)vs, seq_len*sizeof(Point),cudaMemcpyHostToDevice));
// 	cuda_distance<<<seq_len/1024+1, 1024>>>(dist_d, p, vs_d, seq_len,gpu->degree_per_kilometer);
// 	CUDA_SAFE_CALL(cudaMemcpy(&dist, dist_d, sizeof(double), cudaMemcpyDeviceToHost));
// 	gpu->free((char *)dist_d, sizeof(double));
// 	gpu->free((char *)vs_d, seq_len* sizeof(Point));
// 	return dist;
// }



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

