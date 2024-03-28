#include "../include/Ideal.h"

void *rasterization_unit(void *args){
	query_context *ctx = (query_context *)args;
	query_context *gctx = ctx->global_ctx;

	vector<Ideal *> &ideals = *(vector<Ideal *> *)gctx->target;

	// log("thread %d is started",ctx->thread_id);

	while(ctx->next_batch(10)){
		for(int i=ctx->index;i<ctx->index_end;i++){
			struct timeval start = get_cur_time();
			ideals[i]->rasterization(ctx->vpr);
			ctx->report_progress();
		}
	}
	ctx->merge_global();
	return NULL;
}

void process_rasterization(query_context *gctx){

	log("start rasterizing the referred polygons");
	vector<Ideal *> &ideals = *(vector<Ideal *> *)gctx->target;
	assert(ideals.size()>0);
	gctx->index = 0;
	size_t former = gctx->target_num;
	gctx->target_num = ideals.size();

	struct timeval start = get_cur_time();
	pthread_t threads[gctx->num_threads];
	query_context ctx[gctx->num_threads];
	for(int i=0;i<gctx->num_threads;i++){
		ctx[i] = *gctx;
		ctx[i].thread_id = i;
		ctx[i].global_ctx = gctx;
	}

	for(int i=0;i<gctx->num_threads;i++){
		pthread_create(&threads[i], NULL, rasterization_unit, (void *)&ctx[i]);
	}

	for(int i = 0; i < gctx->num_threads; i++ ){
		void *status;
		pthread_join(threads[i], &status);
	}

	//collect partitioning status
	size_t num_partitions = 0;
	size_t num_crosses = 0;
	size_t num_border_partitions = 0;
	size_t num_edges = 0;
	for(auto ideal : ideals){
		num_partitions += ideal->get_num_pixels();
		num_crosses += ideal->get_num_crosses();
		num_border_partitions += ideal->get_num_pixels(BORDER);
		num_edges += ideal->get_num_border_edge();
	}
	logt("IDEALized %d polygons with (%ld)%ld average pixels %.2f average crosses per pixel %.2f edges per pixel", start,
			ideals.size(),
			num_border_partitions/ideals.size(),
			num_partitions/ideals.size(),
			1.0*num_crosses/num_border_partitions,
			1.0*num_edges/num_border_partitions);

	gctx->index = 0;
	gctx->query_count = 0;
	gctx->target_num = former;
}

void preprocess(query_context *gctx){

	vector<Ideal *> target_ideals;
	target_ideals.insert(target_ideals.end(), gctx->source_ideals.begin(), gctx->source_ideals.end());
	target_ideals.insert(target_ideals.end(), gctx->target_ideals.begin(), gctx->target_ideals.end());
	gctx->target = (void *)&target_ideals;
	if(gctx->use_grid){
		process_rasterization(gctx);
	}

	// if(gctx->use_qtree){
	// 	process_qtree(gctx);
	// }

	// if(gctx->use_vector){
	// 	process_convex_hull(gctx);
	// 	process_mer(gctx);
	// 	target_polygons.clear();
	// 	target_polygons.insert(target_polygons.end(), gctx->source_polygons.begin(), gctx->source_polygons.end());
	// 	process_internal_rtree(gctx);
	// }

	// if(gctx->use_geos){
	// 	process_geos(gctx);
	// }

	target_ideals.clear();
	gctx->target = NULL;
}