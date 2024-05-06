#include "cuda_util.h"
#include "../include/Ideal.h"
#include "mygpu.h"

void cuda_create_buffer(query_context *gctx){
	size_t size = 1024 * 1024 * 100;
    size_t size2 = 1024 * 1024 * 200;
    log("CPU momory:");
	gctx->h_status = new uint8_t[size];
    log("\t%.2f MB\tstatus buffer",1.0*size/1024/1024);

	gctx->h_target_status = new uint8_t[size];
    log("\t%.2f MB\ttarget status buffer",1.0*size/1024/1024);

    gctx->h_edges = new pair<double, double>[size];
    log("\t%.2f MB\tedges buffer",1.0*size*sizeof(pair<double, double>)/1024/1024);

    gctx->h_target_edges = new pair<double, double>[size2];
    log("\t%.2f MB\ttarget edges buffer",1.0*size2*sizeof(pair<double, double>)/1024/1024);

	log("GPU memory:");
    CUDA_SAFE_CALL(cudaMalloc((void **) &gctx->d_status, size * sizeof(uint8_t)));
	log("\t%.2f MB\tstatus buffer",1.0*size/1024/1024);

    CUDA_SAFE_CALL(cudaMalloc((void **) &gctx->d_target_status, size * sizeof(uint8_t)));
	log("\t%.2f MB\ttarget status buffer",1.0*size/1024/1024);

    CUDA_SAFE_CALL(cudaMalloc((void **) &gctx->d_edges, 2 * size * sizeof(pair<double, double>)));
	log("\t%.2f MB\tedges buffer",1.0*size*sizeof(pair<double, double>)/1024/1024);
    
    CUDA_SAFE_CALL(cudaMalloc((void **) &gctx->d_target_edges, 2 * size2 * sizeof(pair<double, double>)));
	log("\t%.2f MB\ttarget edges buffer",1.0*size2*sizeof(pair<double, double>)/1024/1024);
}

void preprocess_for_gpu(query_context *gctx){
    size_t size1 = 1024 * 1024 * 100;
    size_t size2 = 1024 * 1024 * 200;
    size_t size = gctx->temp_pair.size() + 1;
    gctx->h_status_offset = new uint[size];
    gctx->h_target_status_offset = new uint[size1];
    gctx->h_edges_offset = new uint[size1];
    gctx->h_target_edges_offset = new uint[size2];

    CUDA_SAFE_CALL(cudaMalloc((void **) &gctx->d_status_offset, size * sizeof(uint)));
    CUDA_SAFE_CALL(cudaMalloc((void **) &gctx->d_target_status_offset, size1 * sizeof(uint)));
    CUDA_SAFE_CALL(cudaMalloc((void **) &gctx->d_edges_offset, size1 * sizeof(uint)));
    CUDA_SAFE_CALL(cudaMalloc((void **) &gctx->d_target_edges_offset, size2 * sizeof(uint)));
    // compact data
    uint pair_id = 0, r_id = 0, t_id = 0, r_point_id = 0, t_point_id = 0;
    for(auto &tp : gctx->temp_pair){
        gctx->h_status_offset[pair_id ++] = r_id;
        Ideal *source = tp.first;
        Ideal *target = tp.second;   
	    vector<int> pxs = source->retrieve_pixels(target->getMBB());
        vector<int> tpxs;
        assert(pxs.size() > 0);
        for(auto px : pxs){
            gctx->h_status[r_id] = source->show_status(px);
            gctx->h_edges_offset[r_id ++] = r_point_id;
            for(int i = 0; i < source->get_num_sequences(px); i ++){
                auto r = source->get_edge_sequence(source->get_offset(px) + i);
                auto pos = r.first;
			    auto len = r.second;
                for(int j = pos; j < len; j ++){
                    gctx->h_edges[r_point_id ++] = make_pair(source->get_boundary()->p[j].x, source->get_boundary()->p[j].y); 
                }
            }
            // h_status每一个代表一个reference polygon pixel
            box bx =  source->get_pixel_box(source->get_x(px), source->get_y(px));
            tpxs = target->retrieve_pixels(&bx);
            assert(tpxs.size() > 0);
            for(auto tpx : tpxs){
                gctx->h_target_status_offset[r_id] = t_id;
                gctx->h_target_status[t_id] = target->show_status(tpx);
                gctx->h_target_edges_offset[t_id ++] = t_point_id;
                for(int i = 0; i < target->get_num_sequences(tpx); i ++){
                    auto r = target->get_edge_sequence(target->get_offset(tpx) + i);
                    auto pos = r.first;
                    auto len = r.second;
                    for(int j = pos; j < len; j ++){
                        gctx->h_target_edges[t_point_id ++] = make_pair(target->get_boundary()->p[j].x, target->get_boundary()->p[j].y); 
                    }
                }
            }
        }

    }


    gctx->h_status_offset[pair_id] = r_id;
    gctx->h_target_status_offset[r_id] = t_id;
    gctx->h_edges_offset[r_id] = r_point_id;
    gctx->h_target_edges_offset[t_id] = t_point_id;

    CUDA_SAFE_CALL(cudaMemcpy(gctx->d_status_offset, gctx->h_status_offset, size * sizeof(uint), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(gctx->d_status, gctx->h_status, size1 * sizeof(uint8_t), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(gctx->d_target_status, gctx->h_target_status, size1 * sizeof(uint8_t), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(gctx->d_edges_offset, gctx->h_edges_offset, size1 * sizeof(uint), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(gctx->d_edges, gctx->h_edges, size1 * sizeof(pair<double, double>), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(gctx->d_target_edges_offset, gctx->h_target_edges_offset, size2 * sizeof(uint), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(gctx->d_target_edges, gctx->h_target_edges, size2 * sizeof(pair<double, double>), cudaMemcpyHostToDevice));
}