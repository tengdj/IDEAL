#include "cuda_util.h"
#include "../include/Ideal.h"
#include "mygpu.h"

void cuda_create_buffer(query_context *gctx){
	size_t size = BUFFER_SIZE;
    log("CPU momory:");
	gctx->h_status = new uint8_t[size];
	memset(gctx->h_status, -1, size * sizeof(uint8_t));
    log("\t%.2f MB\tstatus buffer",1.0*size/1024/1024);

	log("GPU memory:");
    CUDA_SAFE_CALL(cudaMalloc((void **) &gctx->d_status, size * sizeof(uint8_t)));
	log("\t%.2f MB\tstatus buffer",1.0*size/1024/1024);
}

void preprocess_for_gpu(query_context *gctx){
    size_t size = gctx->temp_pair.size() + 1;
    gctx->h_status_offset = new uint[size];
    CUDA_SAFE_CALL(cudaMalloc((void **) &gctx->d_status_offset, size * sizeof(uint)));
    // compact data
    uint pair_id = 0, id = 0;
    for(auto &tp : gctx->temp_pair){
        gctx->h_status_offset[pair_id] = id;

        Ideal *source = tp.first;
        Ideal *target = tp.second;   
	    vector<int> pxs = source->retrieve_pixels(target->getMBB());
        assert(pxs.size() > 0);
        for(auto px : pxs){
            gctx->h_status[id ++] = source->show_status(px);
        }
        pair_id ++;
    }

    gctx->h_status_offset[pair_id] = id;
    CUDA_SAFE_CALL(cudaMemcpy(gctx->d_status_offset, gctx->h_status_offset, size * sizeof(uint), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(gctx->d_status, gctx->h_status, BUFFER_SIZE * sizeof(uint8_t), cudaMemcpyHostToDevice));
}