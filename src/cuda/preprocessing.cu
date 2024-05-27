#include "cuda_util.h"
#include "../include/Ideal.h"
#include "mygpu.h"

void cuda_create_buffer(query_context *gctx){
    cudaSetDevice(1);

	unsigned long long size = BUFFER_SIZE;
    log("CPU momory:");

    gctx->h_info = (Idealinfo*)new char[size / 4ULL];
    log("\t%.2f MB\tideal info buffer",1.0*size/1024/1024/4);
     
	gctx->h_status = new uint8_t[size / 4ULL];
    log("\t%.2f MB\tstatus buffer",1.0*size/1024/1024/4);

    gctx->h_offset = (uint16_t*)new char[size / 2ULL];
    log("\t%.2f MB\toffset buffer",1.0*size/1024/1024/2);

    gctx->h_edge_sequences = (EdgeSeq *)new char[size];
    log("\t%.2f MB\tedge sequences buffer",1.0*size/1024/1024);

    gctx->h_vertices = (Point *)new char[4ULL * size];
    log("\t%.2f MB\tvertices buffer",4.0*size/1024/1024);

	log("GPU memory:");
    CUDA_SAFE_CALL(cudaMalloc((void **) &gctx->d_info, size/4ULL));
	log("\t%.2f MB\tideal info buffer",1.0*size/1024/1024/4);

    CUDA_SAFE_CALL(cudaMalloc((void **) &gctx->d_status, size/4ULL));
	log("\t%.2f MB\tstatus buffer",1.0*size/1024/1024/4);

    CUDA_SAFE_CALL(cudaMalloc((void **) &gctx->d_offset, size/2ULL));
	log("\t%.2f MB\toffset buffer",1.0*size/1024/1024/2);
    
    CUDA_SAFE_CALL(cudaMalloc((void **) &gctx->d_edge_sequences, size));
	log("\t%.2f MB\tedge sequences buffer",1.0*size/1024/1024);

    CUDA_SAFE_CALL(cudaMalloc((void **) &gctx->d_vertices, size * 4ULL));
	log("\t%.2f MB\tvertices buffer",4.0*size/1024/1024);

}

void preprocess_for_gpu(query_context *gctx){
    cudaSetDevice(1);
    // compact data
    uint iidx = 0, sidx = 0, oidx = 0, eidx = 0, vidx = 0;
    for(auto &tp : gctx->ideal_pairs){
        Ideal *source = tp.first;
        int dimx = source->get_dimx(), dimy = source->get_dimy();
        if(source->idealoffset == nullptr){
            source->idealoffset = new IdealOffset{};

            uint info_size = gctx->ideal_pairs.size();
            Idealinfo idealinfo{source->getMBB(), dimx, dimy, source->get_step_x(), source->get_step_y()};
            memcpy(gctx->h_info+iidx, &idealinfo, sizeof(Idealinfo));
            source->idealoffset->info_start = iidx;
            iidx ++;
            source->idealoffset->info_end = iidx;

            uint status_size = (dimx+1)*(dimy+1) / 4 + 1;
            memcpy(gctx->h_status+sidx, source->get_status(), status_size);
            source->idealoffset->status_start = sidx;
            sidx += status_size;
            source->idealoffset->status_end = sidx;

            uint offset_size = (dimx+1)*(dimy+1) + 1;
            memcpy(gctx->h_offset+oidx, source->get_offset(), offset_size * sizeof(uint16_t));
            source->idealoffset->offset_start = oidx;
            oidx += offset_size;
            source->idealoffset->offset_end = oidx;

            uint edge_sequences_size = source->get_len_edge_sequences();
            memcpy(gctx->h_edge_sequences+eidx, source->get_edge_sequence(), edge_sequences_size * sizeof(EdgeSeq));
            source->idealoffset->edge_sequences_start = eidx;
            eidx += edge_sequences_size;
            source->idealoffset->edge_sequences_end = eidx;

            uint vertices_size = source->get_num_vertices();
            memcpy(gctx->h_vertices+vidx, source->get_boundary()->p, vertices_size * sizeof(Point));
            source->idealoffset->vertices_start = vidx;
            vidx += vertices_size;
            source->idealoffset->vertices_end = vidx;
        }

        Ideal *target = tp.second;   
	    dimx = target->get_dimx(), dimy = target->get_dimy();
        if(target->idealoffset == nullptr){
            target->idealoffset = new IdealOffset{};

            uint info_size = gctx->ideal_pairs.size();
            Idealinfo idealinfo{target->getMBB(), dimx, dimy, target->get_step_x(), target->get_step_y()};
            memcpy(gctx->h_info+iidx, &idealinfo, sizeof(Idealinfo));
            target->idealoffset->info_start = iidx;
            iidx ++;
            target->idealoffset->info_end = iidx;

            uint status_size = (dimx+1)*(dimy+1) / 4 + 1;
            assert((status_size+sidx) < 1U * BUFFER_SIZE);
            memcpy(gctx->h_status+sidx, target->get_status(), status_size);
            target->idealoffset->status_start = sidx;
            sidx += status_size;
            target->idealoffset->status_end = sidx;

            uint offset_size = (dimx+1)*(dimy+1) + 1;
            memcpy(gctx->h_offset+oidx, target->get_offset(), offset_size * sizeof(uint16_t));
            target->idealoffset->offset_start = oidx;
            oidx += offset_size;
            target->idealoffset->offset_end = oidx;

            uint edge_sequences_size = target->get_len_edge_sequences();
            assert((edge_sequences_size+eidx)*sizeof(EdgeSeq) < 1U * BUFFER_SIZE);
            memcpy(gctx->h_edge_sequences+eidx, target->get_edge_sequence(), edge_sequences_size * sizeof(EdgeSeq));
            target->idealoffset->edge_sequences_start = eidx;
            eidx += edge_sequences_size;
            target->idealoffset->edge_sequences_end = eidx;

            uint vertices_size = target->get_num_vertices();
            assert((vertices_size+vidx)*sizeof(Point) < (4ULL * BUFFER_SIZE));
            memcpy(gctx->h_vertices+vidx, target->get_boundary()->p, vertices_size * sizeof(Point));
            target->idealoffset->vertices_start = vidx;
            vidx += vertices_size;
            target->idealoffset->vertices_end = vidx;
        }
    }

    CUDA_SAFE_CALL(cudaMemcpy(gctx->d_info, gctx->h_info, BUFFER_SIZE / 4UL * sizeof(uint8_t), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(gctx->d_status, gctx->h_status, BUFFER_SIZE / 4UL * sizeof(uint8_t), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(gctx->d_offset, gctx->h_offset, BUFFER_SIZE / 2UL * sizeof(uint8_t), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(gctx->d_edge_sequences, gctx->h_edge_sequences, BUFFER_SIZE * sizeof(uint8_t), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(gctx->d_vertices, gctx->h_vertices, 4UL * BUFFER_SIZE * sizeof(uint8_t), cudaMemcpyHostToDevice));
}