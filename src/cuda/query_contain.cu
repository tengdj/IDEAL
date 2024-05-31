#include "geometry.cuh"

#define BLOCK_SIZE 1024

struct PixPair{
	int source_pixid = 0;
	int target_pixid = 0;
	int pair_id = 0;
};

struct Batch{
	uint s_start = 0;
	uint t_start = 0;
	uint s_length = 0;
	uint t_length = 0;
	int pair_id = 0;
};

__global__ void kernel_filter(pair<IdealOffset, IdealOffset> *d_pairs, Idealinfo *d_info, uint8_t *d_status, uint size, uint8_t *resultmap, PixPair *d_pixpairs, uint *pp_size){
	const int x = blockIdx.x * blockDim.x + threadIdx.x;
	if(x < size){
		pair<IdealOffset, IdealOffset> &temp_pair = d_pairs[x];
		IdealOffset &source = temp_pair.first;
		IdealOffset &target = temp_pair.second;
		
		const box &s_mbr = d_info[source.info_start].mbr, &t_mbr = d_info[target.info_start].mbr;
		const double &s_step_x = d_info[source.info_start].step_x, &s_step_y = d_info[source.info_start].step_y;
		const int &s_dimx = d_info[source.info_start].dimx, &s_dimy = d_info[source.info_start].dimy;
		const double &t_step_x = d_info[target.info_start].step_x, &t_step_y = d_info[target.info_start].step_y;
		const int &t_dimx = d_info[target.info_start].dimx, &t_dimy = d_info[target.info_start].dimy;

		const int start_x = gpu_get_offset_x(s_mbr.low[0], t_mbr.low[0], s_step_x, s_dimx);
		const int start_y = gpu_get_offset_y(s_mbr.low[1], t_mbr.low[1], s_step_y, s_dimy);
		const int end_x = gpu_get_offset_x(s_mbr.low[0], t_mbr.high[0], s_step_x, s_dimx);
		const int end_y = gpu_get_offset_y(s_mbr.low[1], t_mbr.high[1], s_step_y, s_dimy);

		uint itn = 0, etn = 0;
		for(int i=start_x;i<=end_x;i++){
			for(int j=start_y;j<=end_y;j++){
				int p = gpu_get_id(i , j, s_dimx);
				if(gpu_show_status(d_status, source.status_start, p) == IN) itn ++;
				else if(gpu_show_status(d_status, source.status_start, p) == OUT) etn ++;

				box bx = gpu_get_pixel_box(i, j, s_mbr.low[0], s_mbr.low[1], s_step_x, s_step_y);
				 
				const int _start_x = gpu_get_offset_x(t_mbr.low[0], bx.low[0], t_step_x, t_dimx);
				const int _start_y = gpu_get_offset_y(t_mbr.low[1], bx.low[1], t_step_y, t_dimy);
				const int _end_x = gpu_get_offset_x(t_mbr.low[0], bx.high[0], t_step_x, t_dimx);
				const int _end_y = gpu_get_offset_y(t_mbr.low[1], bx.high[1], t_step_y, t_dimy);

				for(int _i = _start_x; _i <= _end_x; _i ++){
					for(int _j = _start_y; _j <= _end_y; _j ++){
				 		int p2 = gpu_get_id(_i, _j, t_dimx);
						if(gpu_show_status(d_status, source.status_start, p) == OUT && gpu_show_status(d_status, target.status_start, p2) == IN){
						 	resultmap[x] = 2;
							break;
						}
				 		if(gpu_show_status(d_status, source.status_start, p) == OUT && gpu_show_status(d_status, target.status_start, p2) == BORDER){
							int idx = atomicAdd(pp_size, 1U);
							d_pixpairs[idx].source_pixid = p;
							d_pixpairs[idx].target_pixid = p2;
							d_pixpairs[idx].pair_id = x;
				 		
				 		}
						if(gpu_show_status(d_status, source.status_start, p) == BORDER && gpu_show_status(d_status, target.status_start, p2) == BORDER){
							int idx = atomicAdd(pp_size, 1U);
							d_pixpairs[idx].source_pixid = p;
							d_pixpairs[idx].target_pixid = p2;
							d_pixpairs[idx].pair_id = x;
						}
					}
				}
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

__global__ void kernel_ideal_filter(PixPair *d_pixpairs, pair<IdealOffset, IdealOffset> *d_pairs, Idealinfo *d_info, uint8_t *d_status, uint16_t *d_offset, EdgeSeq *d_edge_sequences, Point *d_vertices, uint *size, uint8_t *resultmap){
	const int x = blockIdx.x * blockDim.x + threadIdx.x;
	if(x < *size){
		int p = d_pixpairs[x].source_pixid;
		int p2 = d_pixpairs[x].target_pixid;
		int pair_id = d_pixpairs[x].pair_id;
		if(resultmap[pair_id] != 0) return;		

		pair<IdealOffset, IdealOffset> temp_pair = d_pairs[pair_id];
		IdealOffset source = temp_pair.first;
		IdealOffset target = temp_pair.second;

		const box &s_mbr = d_info[source.info_start].mbr;
		const double &s_step_x = d_info[source.info_start].step_x, &s_step_y = d_info[source.info_start].step_y;
		const int &s_dimx = d_info[source.info_start].dimx, &s_dimy = d_info[source.info_start].dimy;

		if(gpu_show_status(d_status, source.status_start, p) == OUT && gpu_show_status(d_status, target.status_start, p2) == BORDER){
			box bx = gpu_get_pixel_box(gpu_get_x(p, s_dimx), gpu_get_y(p, s_dimx, s_dimy), s_mbr.low[0], s_mbr.low[1], s_step_x, s_step_y);
			Point pix_border[5];
			pix_border[0] = {bx.low[0], bx.low[1]};
			pix_border[1] = {bx.low[0], bx.high[1]};
			pix_border[2] = {bx.high[0], bx.high[1]};
			pix_border[3] = {bx.high[0], bx.low[1]};
			pix_border[4] = {bx.low[0], bx.low[1]};

			const uint &t_offset_start = target.offset_start;
			const uint &t_edge_sequences_start = target.edge_sequences_start; 
			const uint16_t &t_num_sequence = (d_offset+t_offset_start)[p2 + 1] - (d_offset+t_offset_start)[p2];
			const uint &t_vertices_start = target.vertices_start;

			for(int i = 0; i < t_num_sequence; ++ i){
				EdgeSeq r = (d_edge_sequences+t_edge_sequences_start)[(d_offset+t_offset_start)[p2] + i];
				if(segment_intersect_batch(d_vertices+t_vertices_start+r.start, pix_border, r.length, 4)){
					resultmap[pair_id] = 2;
					return;
				}
				
			}
		}
	}
}

__global__ void kernel_unroll(PixPair *d_pixpairs, pair<IdealOffset, IdealOffset> *d_pairs, uint8_t *d_status, uint16_t *d_offset, EdgeSeq *d_edge_sequences, uint *size, Batch *batches, uint *batch_size, uint8_t *resultmap){
	const int x = blockIdx.x * blockDim.x + threadIdx.x;
	if(x < *size){
		int p = d_pixpairs[x].source_pixid;
		int p2 = d_pixpairs[x].target_pixid;
		int pair_id = d_pixpairs[x].pair_id;
		if(resultmap[pair_id] != 0) return;

		pair<IdealOffset, IdealOffset> temp_pair = d_pairs[pair_id];
		IdealOffset source = temp_pair.first;
		IdealOffset target = temp_pair.second;	

		if(gpu_show_status(d_status, source.status_start, p) == BORDER && gpu_show_status(d_status, target.status_start, p2) == BORDER){

			uint s_offset_start = source.offset_start, t_offset_start = target.offset_start; 
			uint s_edge_sequences_start = source.edge_sequences_start, t_edge_sequences_start = target.edge_sequences_start; 
			int s_num_sequence = (d_offset+s_offset_start)[p + 1] - (d_offset+s_offset_start)[p];
			int t_num_sequence = (d_offset+t_offset_start)[p2 + 1] - (d_offset+t_offset_start)[p2];
			uint s_vertices_start = source.vertices_start, t_vertices_start = target.vertices_start;

			for(int i = 0; i < s_num_sequence; ++ i){
				EdgeSeq r = (d_edge_sequences+s_edge_sequences_start)[(d_offset+s_offset_start)[p] + i]; 
				for(int j = 0; j < t_num_sequence; ++ j){
					EdgeSeq r2 = (d_edge_sequences+t_edge_sequences_start)[(d_offset+t_offset_start)[p2] + j];
					int max_size = 32;
					for(uint s = 0; s < r.length; s += max_size){
						uint end_s = min(s + max_size, r.length);
						for(uint t = 0; t < r2.length; t += max_size){
							uint end_t = min(t + max_size, r2.length);
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

		if(segment_intersect_batch((d_vertices+s1), (d_vertices+s2), len1, len2)){
			resultmap[pair_id] = 3;
			return;
		}
	}
}

uint cuda_contain(query_context *gctx){
	cudaSetDevice(1);

    CudaTimer timer;
	
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

	// resultmap status: 0(undecided), 1(contain), 2(not contain)
	uint8_t *d_resultmap = nullptr;
	CUDA_SAFE_CALL(cudaMalloc((void **) &d_resultmap, size * sizeof(uint8_t)));
	CUDA_SAFE_CALL(cudaMemset(d_resultmap, 0, size * sizeof(uint8_t)));

	PixPair *h_pixpairs = new PixPair[32*1024*1024];
	PixPair *d_pixpairs = nullptr;
	CUDA_SAFE_CALL(cudaMalloc((void **) &d_pixpairs, 32*1024*1024*sizeof(PixPair)));
	
	uint *d_pp_size = nullptr;
	CUDA_SAFE_CALL(cudaMalloc((void **) &d_pp_size, sizeof(uint)));
	CUDA_SAFE_CALL(cudaMemset(d_pp_size, 0, sizeof(uint)));

	/*1. Raster Model Filtering*/

	int grid_size_x = (size + BLOCK_SIZE - 1) / BLOCK_SIZE;
	dim3 block_size(BLOCK_SIZE, 1, 1);
	dim3 grid_size(grid_size_x, 1, 1);
	
    timer.startTimer();

	kernel_filter<<<grid_size, block_size>>>(d_pairs, gctx->d_info, gctx->d_status, size, d_resultmap, d_pixpairs, d_pp_size);
	cudaDeviceSynchronize();
	check_execution();
	
	timer.stopTimer();
    printf("kernel_filter time: %f ms\n", timer.getElapsedTime());

	uint h_pp_size;
	CUDA_SAFE_CALL(cudaMemcpy(&h_pp_size, d_pp_size, sizeof(uint), cudaMemcpyDeviceToHost));

	printf("h_pp_size: %u\n", h_pp_size);

	// int fd = 0;
	// for(int i = 0; i < size; ++ i ){
	// 	if(h_resultmap[i] == 1 || h_resultmap[i] == 0) fd ++;
	// }
	// return fd;
	return 0;

	/*2. Ideal Filtering */

	grid_size_x = (h_pp_size + BLOCK_SIZE - 1) / BLOCK_SIZE;
	grid_size.x = grid_size_x;

	timer.startTimer();

	kernel_ideal_filter<<<grid_size, block_size>>>(d_pixpairs, d_pairs, gctx->d_info, gctx->d_status, gctx->d_offset, gctx->d_edge_sequences, gctx->d_vertices, d_pp_size, d_resultmap);
	cudaDeviceSynchronize();
    check_execution();

	timer.stopTimer();
    printf("Kernel_ideal_filter execution time: %f ms\n", timer.getElapsedTime());

	/*3. Unroll Refinement*/

	Batch *d_batch = nullptr;
	CUDA_SAFE_CALL(cudaMalloc((void **) &d_batch, 64 * 1024 * 1024 * sizeof(Batch)));
	uint *d_batch_size = nullptr;
	CUDA_SAFE_CALL(cudaMalloc((void **) &d_batch_size, sizeof(uint)));
	CUDA_SAFE_CALL(cudaMemset(d_batch_size, 0, sizeof(uint)));
	
    timer.startTimer();

	kernel_unroll<<<grid_size, block_size>>>(d_pixpairs, d_pairs, gctx->d_status, gctx->d_offset, gctx->d_edge_sequences, d_pp_size, d_batch, d_batch_size, d_resultmap);
	cudaDeviceSynchronize();
    check_execution();

	timer.stopTimer();
    printf("Kernel_unroll execution time: %f ms\n", timer.getElapsedTime());


	uint h_batch_size;
	CUDA_SAFE_CALL(cudaMemcpy(&h_batch_size, d_batch_size, sizeof(uint), cudaMemcpyDeviceToHost));

	// cout << "Batch Size = " << h_batch_size << endl;

	/*4. Refinement step*/

	grid_size_x = (h_batch_size + BLOCK_SIZE - 1) / BLOCK_SIZE;
	grid_size.x = grid_size_x;

	timer.startTimer();

	kernel_refinement<<<grid_size, block_size>>>(d_batch, gctx->d_vertices, d_batch_size, d_resultmap);
	cudaDeviceSynchronize();
	check_execution();

	timer.stopTimer();
    printf("Kernel_refinement execution time: %f ms\n", timer.getElapsedTime());

	uint8_t *h_resultmap = new uint8_t[size];
	CUDA_SAFE_CALL(cudaMemcpy(h_resultmap, d_resultmap, size * sizeof(uint8_t), cudaMemcpyDeviceToHost));

	int found = 0;
	for(int i = 0; i < size; ++ i ){
		if(h_resultmap[i] == 1 || h_resultmap[i] == 0) found ++;
	}

	delete []h_pairs;
	delete []h_resultmap;
	delete []h_pixpairs;

	CUDA_SAFE_CALL(cudaFree(d_pairs));
	CUDA_SAFE_CALL(cudaFree(d_resultmap));
	CUDA_SAFE_CALL(cudaFree(d_pixpairs));
	CUDA_SAFE_CALL(cudaFree(d_batch));
	CUDA_SAFE_CALL(cudaFree(d_batch_size));

	return found;
}
