#include "cuda_util.h"
#include "Ideal.h"

__device__ int gpu_get_id(int x, int y, int dimx){
	return y * (dimx+1) + x;
}

// from id to pixel x
__device__ int gpu_get_x(int id, int dimx){
	return id % (dimx+1);
}

// from id to pixel y
__device__ int gpu_get_y(int id, int dimx, int dimy){
	assert((id / (dimx+1)) <= dimy);
	return id / (dimx+1);
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

__device__ PartitionStatus gpu_show_status(uint8_t *status, uint &start, int &id){
	uint8_t st = (status+start)[id / 4];
	int pos = id % 4 * 2;   // The multiplication by 2 is because each status occupies 2 bits.	
	st &= ((uint8_t)3 << pos);
	st >>= pos;
	if(st == 0) return OUT;
	if(st == 3) return IN;
	return BORDER;
}

__device__ box gpu_get_pixel_box(int x, int y, double bx_lowx, double bx_lowy, double step_x, double step_y){
	const double start_x = bx_lowx;
	const double start_y = bx_lowy;

	double lowx = start_x + x * step_x;
	double lowy = start_y + y * step_y;
	double highx = start_x + (x + 1) * step_x;
	double highy = start_y + (y + 1) * step_y;

	return box(lowx, lowy, highx, highy);
}

// __device__ bool gpu_segment_intersect(const Point& a, const Point& b, box &bx){
// 	if (c.cross(a, d) == 0 && c.cross(b, d) == 0)
//         return inter1(a.x, b.x, c.x, d.x) && inter1(a.y, b.y, c.y, d.y);
//     return sgn(a.cross(b, c)) != sgn(a.cross(b, d)) &&
//            sgn(c.cross(d, a)) != sgn(c.cross(d, b));
// }

// __device__ bool gpu_segment_intersect_box(Point* p, uint s, box &bx){
// 	for(int i = 0; i < s; i ++){
// 		if(gpu_segment_intersect(p[i], p[i+1], bx)) 
// 			return true;
// 	}
// 	return false;
// }