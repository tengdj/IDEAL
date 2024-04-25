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
    printf("temp_pair.size()=%d\n", gctx->temp_pair.size());
    // compact data
    uint pair_id = 0, id = 0;
    unsigned long long test_id = 0;
    for(auto &tp : gctx->temp_pair){
        Ideal *source = tp.first;
        Ideal *target = tp.second;
        if(!source->getMBB()->contain(*target->getMBB())){
            continue;
        }
        
        gctx->h_status_offset[pair_id] = id;
	    vector<int> pxs = source->retrieve_pixels(target->getMBB());
        test_id += pxs.size();
        assert(pxs.size() > 0);
        for(auto px : pxs){
            if(source->show_status(px) == IN) gctx->h_status[id ++] = 0;
            if(source->show_status(px) == BORDER) gctx->h_status[id ++] = 1;
            if(source->show_status(px) == OUT) gctx->h_status[id ++] = 2;
        }
        pair_id ++;
    }
    printf("TEST ID: %d\n", test_id);

    gctx->h_status_offset[pair_id] = id;
    printf("pair_id = %d, id=%d\n", pair_id, id);
    
    CUDA_SAFE_CALL(cudaMemcpy(gctx->d_status_offset, gctx->h_status_offset, size * sizeof(uint), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(gctx->d_status, gctx->h_status, BUFFER_SIZE * sizeof(uint8_t), cudaMemcpyHostToDevice));
}