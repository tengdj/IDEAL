#include "cuda_util.h"
#include "Ideal.h"
#include "mygpu.h"

#define BLOCK_SIZE 1024

struct Threadwork{
	int source_pixid = 0;
	int target_pixid = 0;
	int pair_id = 0;
};

struct Batch{
	uint s_start = 0;
	uint t_start = 0;
	uint8_t s_length = 0;
	uint8_t t_length = 0;
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

__device__ int gpu_get_offset_x(double s_xval, double t_xval, double step_x, int dimx){
	int x = gpu_double_to_int((t_xval-s_xval)/step_x);
	return min(max(x, 0), dimx);
}

__device__ int gpu_get_offset_y(double s_yval, double t_yval, double step_y, int dimy){
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


__global__ void kernel_filter(pair<IdealOffset, IdealOffset> *d_pairs,Idealinfo *d_info, uint8_t *d_status, uint size, uint8_t *resultmap){
	const int x = blockIdx.x * blockDim.x + threadIdx.x;
	if(x < size){
		pair<IdealOffset, IdealOffset> temp_pair = d_pairs[x];
		IdealOffset source = temp_pair.first;
		IdealOffset target = temp_pair.second;		
		
		box s_mbr = d_info[source.info_start].mbr, t_mbr = d_info[target.info_start].mbr;
		double s_step_x = d_info[source.info_start].step_x, s_step_y = d_info[source.info_start].step_y;
		int s_dimx = d_info[source.info_start].dimx, s_dimy = d_info[source.info_start].dimy;

		int start_x = gpu_get_offset_x(s_mbr.low[0], t_mbr.low[0], s_step_x, s_dimx);
		int start_y = gpu_get_offset_y(s_mbr.low[1], t_mbr.low[1], s_step_y, s_dimy);
		int end_x = gpu_get_offset_x(s_mbr.low[0], t_mbr.high[0], s_step_x, s_dimx);
		int end_y = gpu_get_offset_y(s_mbr.low[1], t_mbr.high[1], s_step_y, s_dimy);

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

__global__ void process_heavy_workload(Point *p1, Point *p2, uint s1, uint s2, int pair_id, uint8_t* resultmap){
	const int x = blockIdx.x * blockDim.x + threadIdx.x;
    const int y = blockIdx.y * blockDim.y + threadIdx.y;
	if(x < s1 && y < s2){
		if(resultmap[pair_id] != 0) return;
		if(gpu_segment_intersect(p1[x],p1[(x+1)%s1],p2[y],p2[(y+1)%s2])){
			resultmap[pair_id] = 3;
			return;
		}
	}
}

__global__ void kernel_unroll(Threadwork *d_threadwork, pair<IdealOffset, IdealOffset> *d_pairs, uint16_t *d_offset, EdgeSeq *d_edge_sequences, Point *d_vertices, uint size, Batch *batches, uint *batch_size, uint8_t *resultmap){
	const int x = blockIdx.x * blockDim.x + threadIdx.x;
	if(x < size){
		int p = d_threadwork[x].source_pixid;
		int p2 = d_threadwork[x].target_pixid;
		int pair_id = d_threadwork[x].pair_id;
		// if(resultmap[pair_id] != 0) return;

		pair<IdealOffset, IdealOffset> temp_pair = d_pairs[pair_id];
		IdealOffset source = temp_pair.first;
		IdealOffset target = temp_pair.second;	

		uint s_offset_start = source.offset_start, t_offset_start = target.offset_start; 
		uint s_edge_sequences_start = source.edge_sequences_start, t_edge_sequences_start = target.edge_sequences_start; 
		int s_num_sequence = (d_offset+s_offset_start)[p + 1] - (d_offset+s_offset_start)[p];
		int t_num_sequence = (d_offset+t_offset_start)[p2 + 1] - (d_offset+t_offset_start)[p2];
		uint s_vertices_start = source.vertices_start, t_vertices_start = target.vertices_start;
		// int idx = atomicAdd(batch_size, 1U*s_num_sequence*t_num_sequence);
		for(int i = 0; i < s_num_sequence; ++ i){
			EdgeSeq r = (d_edge_sequences+s_edge_sequences_start)[(d_offset+s_offset_start)[p] + i]; 
			for(int j = 0; j < t_num_sequence; ++ j){
				EdgeSeq r2 = (d_edge_sequences+t_edge_sequences_start)[(d_offset+t_offset_start)[p2] + j];
				if(r.length <= 32 && r2.length <= 32){
					if(gpu_segment_intersect_batch((d_vertices+s_vertices_start+r.start), (d_vertices+t_vertices_start+r2.start), r.length, r2.length)){
						resultmap[pair_id] = 3;
						return;
					}
				}else{
					for(int s = 0; s < r.length; s += 32){
						for(int t = 0; t < r2.length; t += 32){
							int end_s = min(s + 32, r.length);
							int end_t = min(t + 32, r2.length);
							uint idx = atomicAdd(batch_size, 1U);
							batches[idx].s_start = s_vertices_start+r.start + s;
							batches[idx].t_start = t_vertices_start+r2.start + t;
							batches[idx].s_length = end_s - s;
							batches[idx].t_length = end_t - t;
							batches[idx].pair_id = pair_id;
						}
					}	
				}
			}
		}
	}
}

__global__ void kernel_refinement(Batch* batches, Point *d_vertices, uint *size, uint8_t* resultmap){
	const int x = blockIdx.x * blockDim.x + threadIdx.x;
	if(x < *size){
		uint s1  = batches[x].s_start;
		uint s2 = batches[x].t_start;
		uint len1 = batches[x].s_length;
		uint len2 = batches[x].t_length;
		int pair_id = batches[x].pair_id;
		if(resultmap[pair_id] != 0) return;

		// if(len1 * len2 > 10000) {
		// 	const int grid_size_x = (len1 + 32 - 1) / 32;
    	// 	const int grid_size_y = (len2 + 32 - 1) / 32;
		// 	const dim3 block_size(32, 32);
    	// 	const dim3 grid_size(grid_size_x, grid_size_y);
		// 	process_heavy_workload<<<grid_size, block_size>>>(d_vertices+s1, d_vertices+s2, len1, len2, pair_id, resultmap);
		// 	cudaDeviceSynchronize();
		// }
		if(gpu_segment_intersect_batch((d_vertices+s1), (d_vertices+s2), len1, len2)){
			resultmap[pair_id] = 3;
			return;
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
	// 定义事件变量
    cudaEvent_t start, stop;
    float elapsedTime;
	
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
	
	 // 创建事件
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
	// 记录起始事件
    cudaEventRecord(start, 0);
	kernel_filter<<<grid_size, block_size>>>(d_pairs, gctx->d_info, gctx->d_status, size, d_resultmap);
	cudaDeviceSynchronize();
	cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess) {
        printf("FILTER CUDA error: %s\n", cudaGetErrorString(error));
    }
	// 记录结束事件
    cudaEventRecord(stop, 0);
	// 同步事件
    cudaEventSynchronize(stop);
	// 计算时间差
    cudaEventElapsedTime(&elapsedTime, start, stop);

    // 打印运行时间
    printf("kernel_filter time: %f ms\n", elapsedTime);
	

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

	Batch *d_batch = nullptr;
	CUDA_SAFE_CALL(cudaMalloc((void **) &d_batch, 1024 * 1024 * 1024 * sizeof(Batch)));
	uint *d_batch_size = nullptr;
	CUDA_SAFE_CALL(cudaMalloc((void **) &d_batch_size, sizeof(uint)));
	CUDA_SAFE_CALL(cudaMemset(d_batch_size, 0, sizeof(uint)));

	int id = atomic_id.load();

	CUDA_SAFE_CALL(cudaMemcpy(d_threadwork, h_threadwork, 8*1024*1024*sizeof(Threadwork), cudaMemcpyHostToDevice));

	grid_size_x = (id + BLOCK_SIZE - 1) / BLOCK_SIZE;
	grid_size.x = grid_size_x;
	
    // 创建事件
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
	// 记录起始事件
    cudaEventRecord(start, 0);
	kernel_unroll<<<grid_size, block_size>>>(d_threadwork, d_pairs, gctx->d_offset, gctx->d_edge_sequences, gctx->d_vertices, id, d_batch, d_batch_size, d_resultmap);
	// 记录结束事件
    cudaEventRecord(stop, 0);
	// 同步事件
    cudaEventSynchronize(stop);
	// 计算时间差
    cudaEventElapsedTime(&elapsedTime, start, stop);

    // 打印运行时间
    printf("Kernel_unroll execution time: %f ms\n", elapsedTime);


	uint h_batch_size;
	CUDA_SAFE_CALL(cudaMemcpy(&h_batch_size, d_batch_size, sizeof(uint), cudaMemcpyDeviceToHost));

	cout << "Batch Size = " << h_batch_size << endl;

	grid_size_x = (h_batch_size + BLOCK_SIZE - 1) / BLOCK_SIZE;
	grid_size.x = grid_size_x;
	// 记录起始事件
    cudaEventRecord(start, 0);
	kernel_refinement<<<grid_size, block_size>>>(d_batch, gctx->d_vertices, d_batch_size, d_resultmap);
	cudaDeviceSynchronize();
	// 记录结束事件
    cudaEventRecord(stop, 0);
	// 同步事件
    cudaEventSynchronize(stop);
	// 计算时间差
    cudaEventElapsedTime(&elapsedTime, start, stop);

    // 打印运行时间
    printf("Kernel_refinement execution time: %f ms\n", elapsedTime);

	error = cudaGetLastError();
    if (error != cudaSuccess) {
        printf("REFINEMENT CUDA error: %s\n", cudaGetErrorString(error));
    }

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

