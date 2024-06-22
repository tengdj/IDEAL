#include "geometry.cuh"

#define WITHIN_DISTANCE 10

struct PixPolyPair{
    
};

__global__ void kernel_init(pair<IdealOffset, IdealOffset> *d_pairs, Idealinfo *d_info, uint size, double *distance){
	const int x = blockIdx.x * blockDim.x + threadIdx.x;
	if(x < size){
		pair<IdealOffset, IdealOffset> &pair = d_pairs[x];
		IdealOffset &source = pair.first;
		IdealOffset &target = pair.second;
		box &s_mbr = d_info[source.info_start].mbr;
        box &t_mbr = d_info[target.info_start].mbr;

		distance[x] = gpu_max_distance(s_mbr, t_mbr);
	}
}

__global__ void kernel_1(pair<IdealOffset, IdealOffset> *d_pairs, Idealinfo *d_info, uint8_t *d_status, PixMapping* d_pixpolypairs, uint *buffer_size, uint size, int *step, bool *resultmap){
    const int pair_id = blockIdx.x * blockDim.x + threadIdx.x;
    if(pair_id < size){
		if(resultmap[pair_id] != 0) return;
        pair<IdealOffset, IdealOffset> &pair = d_pairs[pair_id];
        IdealOffset &source = pair.first;
        IdealOffset &target = pair.second;

        box &s_mbr = d_info[source.info_start].mbr, &t_mbr = d_info[target.info_start].mbr;
        const double &s_step_x = d_info[source.info_start].step_x, &s_step_y = d_info[source.info_start].step_y;
		const int &s_dimx = d_info[source.info_start].dimx, &s_dimy = d_info[source.info_start].dimy;

        int lowx = gpu_get_offset_x(s_mbr.low[0], t_mbr.low[0], s_step_x, s_dimx);
        int lowy = gpu_get_offset_x(s_mbr.low[1], t_mbr.low[1], s_step_y, s_dimy);
        int highx = gpu_get_offset_y(s_mbr.low[0], t_mbr.high[0], s_step_x, s_dimx);
        int highy = gpu_get_offset_y(s_mbr.low[1], t_mbr.high[1], s_step_y, s_dimy);

        int ymin = max(0, lowy - *step);
        int ymax = min(s_dimy, highy + *step);

        // left scan      
        if(lowx - *step >= 0){
            for(int y = ymin; y <= ymax; y ++){
                int id = gpu_get_id(lowx-*step, y, s_dimx);
                if(gpu_show_status(d_status, source.status_start, id) == BORDER){
                    int idx = atomicAdd(buffer_size, 1);
                    d_pixpolypairs[idx] = {pair_id, id};
                }
            }
        }
        // right scan
        if(highx + *step <= s_dimx){
            for(int y = ymin; y <= ymax; y++){
                int id = gpu_get_id(highx+*step, y, s_dimx);
                if(gpu_show_status(d_status, source.status_start, id) == BORDER){
                    int idx = atomicAdd(buffer_size, 1);
                    d_pixpolypairs[idx] = {pair_id, id};
                }
            }
        }
        // skip the first if there is left scan
        int xmin = max(0, lowx - *step + (lowx - *step >= 0));
        // skip the last if there is right scan
        int xmax = min(s_dimx, highx + *step - (highx + *step <= s_dimx));  

        // bottom scan
        if (lowy - *step >= 0) {
            for (int x = xmin; x <= xmax; x++) {
                int id = gpu_get_id(x, lowy-*step, s_dimx);
                if(gpu_show_status(d_status, source.status_start, id) == BORDER){
                    int idx = atomicAdd(buffer_size, 1);
                    d_pixpolypairs[idx] = {pair_id, id};
                }
            }
        }   
        // top scan
        if (highy + *step <= s_dimy) {
            for (int x = xmin; x <= xmax; x++) {
                int id = gpu_get_id(x, highy+*step, s_dimx);
                if(gpu_show_status(d_status, source.status_start, id) == BORDER){
                    int idx = atomicAdd(buffer_size, 1);
                    d_pixpolypairs[idx] = {pair_id, id};
                }
            }
        }
    }    
}

__global__ void kernel_2_1(PixMapping *d_pixpolypairs, pair<IdealOffset, IdealOffset> *d_pairs, Idealinfo *d_info, uint8_t *d_status, int *step, uint *size, PixPair *d_pixpairs, uint *buffer_size, bool *resultmap){
    const int x = blockIdx.x * blockDim.x + threadIdx.x;
    if(x < *size){
		int pair_id = d_pixpolypairs[x].pair_id;
		int cur = d_pixpolypairs[x].pix_id;
        if(resultmap[pair_id]) return;

		pair<IdealOffset, IdealOffset> &pair = d_pairs[pair_id];
        IdealOffset &source = pair.first;
        IdealOffset &target = pair.second;

        box &s_mbr = d_info[source.info_start].mbr, &t_mbr = d_info[target.info_start].mbr;
		const double &s_step_x = d_info[source.info_start].step_x, &s_step_y = d_info[source.info_start].step_y;
		const int &s_dimx = d_info[source.info_start].dimx, &s_dimy = d_info[source.info_start].dimy;
        const double &t_step_x = d_info[target.info_start].step_x, &t_step_y = d_info[target.info_start].step_x;
        const int &t_dimx = d_info[target.info_start].dimx, &t_dimy = d_info[target.info_start].dimy;

        auto pix_box = gpu_get_pixel_box(gpu_get_x(cur, s_dimx), gpu_get_y(cur, s_dimx, s_dimy), s_mbr.low[0], s_mbr.low[1], s_step_x, s_step_y);

        int lowx = gpu_get_offset_x(t_mbr.low[0], pix_box.low[0], t_step_x, t_dimx);
        int lowy = gpu_get_offset_x(t_mbr.low[0], pix_box.high[0], t_step_x, t_dimx);
        int highx = gpu_get_offset_y(t_mbr.low[1], pix_box.low[1], t_step_y, t_dimy);
        int highy = gpu_get_offset_y(t_mbr.low[1], pix_box.high[1], t_step_y, t_dimy);

        int ymin = max(0, lowy - *step);
        int ymax = min(t_dimy, highy + *step);

        // left scan      
        if(lowx - *step >= 0){
            for(int y = ymin; y <= ymax; y ++){
                int id = gpu_get_id(lowx-*step, y, t_dimx);
                if(gpu_show_status(d_status, target.status_start, id) == BORDER){
                    int idx = atomicAdd(buffer_size, 1);
                    d_pixpairs[idx] = {cur, id, pair_id};
                }
            }
        }
        // right scan
        if(highx + *step <= t_dimx){
            for(int y = ymin; y <= ymax; y++){
                int id = gpu_get_id(highx+*step, y, t_dimx);
                if(gpu_show_status(d_status, target.status_start, id) == BORDER){
                    int idx = atomicAdd(buffer_size, 1);
                    d_pixpairs[idx] = {cur, id, pair_id};
                }
            }
        }
        // skip the first if there is left scan
        int xmin = max(0, lowx - *step + (lowx - *step >= 0));
        // skip the last if there is right scan
        int xmax = min(t_dimx, highx + *step - (highx + *step <= t_dimx));  

        // bottom scan
        if (lowy - *step >= 0) {
            for (int x = xmin; x <= xmax; x++) {
                int id = gpu_get_id(x, lowy-*step, t_dimx);
                if(gpu_show_status(d_status, target.status_start, id) == BORDER){
                    int idx = atomicAdd(buffer_size, 1);
                    d_pixpairs[idx] = {cur, id, pair_id};
                }
            }
        }   
        // top scan
        if (highy + *step <= t_dimy) {
            for (int x = xmin; x <= xmax; x++) {
                int id = gpu_get_id(x, highy+*step, t_dimx);
                if(gpu_show_status(d_status, target.status_start, id) == BORDER){
                    int idx = atomicAdd(buffer_size, 1);
                    d_pixpairs[idx] = {cur, id, pair_id};
                }
            }
        }

	}	
}

__global__ void kernel_2_2(PixPair *d_pixpairs, pair<Idealinfo, Idealinfo> *d_pairs, uint16_t *d_offset, EdgeSeq *d_edge_sequences, Point *d_vertices, uint *size, double *distance, bool *resultmap){
    const int x = blockIdx.x * blockDim.x + threadIdx.x;
    if(x < *size){
        int pair_id = d_pixpairs[x].pair_id;
    }
}

uint cuda_within_polygon(query_context *gctx){
    CudaTimer timer;

    float sum_filter = 0.0;
	float sum_unroll = 0.0;
	float sum_refinement = 0.0;
	float sum_check = 0.0;

    uint size = gctx->polygon_pairs.size();

    pair<IdealOffset, IdealOffset> *h_pairs = new pair<IdealOffset, IdealOffset>[size];
    pair<IdealOffset, IdealOffset> *d_pairs = nullptr;

    for(int i = 0; i < size; i ++ ){
        Ideal *source = gctx->polygon_pairs[i].first;
        Ideal *target = gctx->polygon_pairs[i].second;
        h_pairs[i] = {*source->idealoffset, *target->idealoffset};
    }

    CUDA_SAFE_CALL(cudaMalloc((void**) &d_pairs, size * sizeof(pair<IdealOffset, IdealOffset>)));
	CUDA_SAFE_CALL(cudaMemcpy(d_pairs, h_pairs, size * sizeof(pair<IdealOffset, IdealOffset>), cudaMemcpyHostToDevice));
    
	double *h_distance = new double[size * sizeof(double)];
    double *d_distance = nullptr;
    CUDA_SAFE_CALL(cudaMalloc((void **) &d_distance, size * sizeof(double)));

	bool *d_resultmap = nullptr;
    CUDA_SAFE_CALL(cudaMalloc((void **) &d_resultmap, size * sizeof(bool)));
	CUDA_SAFE_CALL(cudaMemset(d_resultmap, 0, size * sizeof(bool)));

    int h_step = 0;
	int *d_step = nullptr;
	CUDA_SAFE_CALL(cudaMalloc((void **) &d_step, sizeof(int)));
	CUDA_SAFE_CALL(cudaMemset(d_step, 0, sizeof(int)));

    int h_step_inner = 0;
    int *d_step_inner = nullptr;
    CUDA_SAFE_CALL(cudaMalloc((void **) &d_step_inner, sizeof(int)));
	CUDA_SAFE_CALL(cudaMemset(d_step_inner, 0, sizeof(int)));

    char *d_BufferInput = nullptr;
    CUDA_SAFE_CALL(cudaMalloc((void **) &d_BufferInput, 4UL * 1024 * 1024 * 1024));
    uint *d_bufferinput_size = nullptr;
    CUDA_SAFE_CALL(cudaMalloc((void **) &d_bufferinput_size, sizeof(uint)));
    uint h_bufferinput_size;

    char *d_BufferOutput = nullptr;
    CUDA_SAFE_CALL(cudaMalloc((void **) &d_BufferOutput, 4UL * 1024 * 1024 * 1024));
    uint *d_bufferoutput_size = nullptr;
    CUDA_SAFE_CALL(cudaMalloc((void **) &d_bufferoutput_size, sizeof(uint)));
    uint h_bufferoutput_size;

    int grid_size_x = (size + BLOCK_SIZE - 1) / BLOCK_SIZE;
	dim3 block_size(BLOCK_SIZE, 1, 1);
	dim3 grid_size(grid_size_x, 1, 1);

    timer.startTimer();

	kernel_init<<<grid_size, block_size>>>(d_pairs, gctx->d_info, size, d_distance);
	cudaDeviceSynchronize();
	check_execution("kernel init");

    timer.stopTimer();
    printf("kernel initialization time: %f ms\n", timer.getElapsedTime());

    while(true){
        printf("STEP: %d\n", h_step);
		CUDA_SAFE_CALL(cudaMemset(d_bufferinput_size, 0, sizeof(uint)));

        timer.startTimer();

		kernel_1<<<grid_size, block_size>>>(d_pairs, gctx->d_info, gctx->d_status, (PixMapping*)d_BufferInput, d_bufferinput_size, size, d_step, d_resultmap);
		cudaDeviceSynchronize();
		check_execution("Kernel filter");

        timer.stopTimer();
    	printf("kernel_1 time: %f ms\n", timer.getElapsedTime());
		sum_filter += timer.getElapsedTime();

		CUDA_SAFE_CALL(cudaMemcpy(&h_bufferinput_size, d_bufferinput_size, sizeof(uint), cudaMemcpyDeviceToHost));

        if(h_bufferinput_size == 0) break;	
        printf("h_bufferinput_size: %u\n", h_bufferinput_size);
        if(h_step == 100) break;
        h_step ++;
        CUDA_SAFE_CALL(cudaMemcpy(d_step, &h_step, sizeof(int), cudaMemcpyHostToDevice));

        while(true){
            printf("STEP_INNER: %d\n", h_step_inner);
		    CUDA_SAFE_CALL(cudaMemset(d_bufferoutput_size, 0, sizeof(uint)));


            grid_size.x = (h_bufferinput_size + BLOCK_SIZE - 1) / BLOCK_SIZE;

            kernel_2_1<<<grid_size, block_size>>>((PixMapping*)d_BufferInput, d_pairs, gctx->d_info, gctx->d_status, d_step_inner, d_bufferinput_size, (PixPair*)d_BufferOutput, d_bufferoutput_size, d_resultmap);
            cudaDeviceSynchronize();
		    check_execution("Kernel_2_1");

            
            CUDA_SAFE_CALL(cudaMemcpy(&h_bufferoutput_size, d_bufferoutput_size, sizeof(uint), cudaMemcpyDeviceToHost));
            printf("h_bufferoutput_size: %u\n", h_bufferoutput_size);

            grid_size.x = (h_bufferoutput_size + BLOCK_SIZE - 1) / BLOCK_SIZE;

            kernel_2_2<<<grid_size, block_size>>>((PixPair*)d_BufferOutput, d_pairs, gctx->d_info, gctx->d_offset, gctx->d_edge_sequences, gctx->d_vertices, d_bufferoutput_size, d_distance, d_resultmap);
            cudaDeviceSynchronize();
		    check_execution("Kernel_2_1");
           
            if(h_step_inner == 10) break;
            h_step_inner ++;
            CUDA_SAFE_CALL(cudaMemcpy(d_step_inner, &h_step_inner, sizeof(int), cudaMemcpyHostToDevice));

            




        }
        h_step_inner = 0;
        CUDA_SAFE_CALL(cudaMemset(d_step_inner, 0, sizeof(int)));

    }

}