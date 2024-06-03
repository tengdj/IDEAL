#pragma once

#include "cuda_util.h"
#include "Ideal.h"

#define BLOCK_SIZE 1024

struct PointPixPair{
	int pair_id = 0;
	int pix_id = 0;
};

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

__device__ __forceinline__ int gpu_get_id(int x, int y, int dimx){
	return y * (dimx+1) + x;
}

// from id to pixel x
__device__ __forceinline__ int gpu_get_x(int id, int dimx){
	return id % (dimx+1);
}

// from id to pixel y
__device__ __forceinline__ int gpu_get_y(int id, int dimx, int dimy){
	assert((id / (dimx+1)) <= dimy);
	return id / (dimx+1);
}

__device__ __forceinline__ int  gpu_double_to_int(double val){
	int vi = (int)val;
	if(abs(1.0*(vi+1)-val)<0.00000001){
		vi++;
	}
	return vi;
}

__device__ __forceinline__ int gpu_get_offset_x(double s_xval, double t_xval, double step_x, int dimx){
	int x = gpu_double_to_int((t_xval-s_xval)/step_x);
	return min(max(x, 0), dimx);
}

__device__ __forceinline__ int gpu_get_offset_y(double s_yval, double t_yval, double step_y, int dimy){
	int y = gpu_double_to_int((t_yval-s_yval)/step_y);
	return min(max(y, 0), dimy);
}

__device__ __forceinline__ PartitionStatus gpu_show_status(uint8_t *status, uint &start, int &id){
	uint8_t st = (status+start)[id / 4];
	int pos = id % 4 * 2;   // The multiplication by 2 is because each status occupies 2 bits.	
	st &= ((uint8_t)3 << pos);
	st >>= pos;
	if(st == 0) return OUT;
	if(st == 3) return IN;
	return BORDER;
}

__device__ __forceinline__ box gpu_get_pixel_box(int x, int y, double bx_lowx, double bx_lowy, double step_x, double step_y){
	const double start_x = bx_lowx;
	const double start_y = bx_lowy;

	double lowx = start_x + x * step_x;
	double lowy = start_y + y * step_y;
	double highx = start_x + (x + 1) * step_x;
	double highy = start_y + (y + 1) * step_y;

	return box(lowx, lowy, highx, highy);
}

__device__ __forceinline__ int gpu_get_pixel_id(Point &p, box &s_mbr, double step_x, double step_y, int dimx, int dimy){
	int xoff = gpu_get_offset_x(s_mbr.low[0], p.x, step_x, dimx);
	int yoff = gpu_get_offset_y(s_mbr.low[1], p.y, step_y, dimy);
	assert(xoff <= dimx);
	assert(yoff <= dimy);
	return gpu_get_id(xoff, yoff, dimx);
}