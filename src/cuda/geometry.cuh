#pragma once

#include "cuda_util.h"
#include "Ideal.h"

#define BLOCK_SIZE 1024

const double EARTH_RADIUS_KM = 6371.0;

struct PointPixPair{
	int pair_id = 0;
	int pix_id = 0;
};

struct PixPair{
	int source_pixid = 0;
	int target_pixid = 0;
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

__device__ __forceinline__ int gpu_get_closest_pixel(Point &p, int pixx, int pixy, int dimx, int dimy){
	if(pixx < 0){
		pixx = 0;
	}
	if(pixx > dimx){
		pixx = dimx;
	}
	if(pixy < 0){
		pixy = 0;
	}
	if(pixy > dimy){
		pixy = dimy;
	}
	return gpu_get_id(pixx, pixy, dimx);	
}

__device__ __forceinline__ double degreesToRadians(double degrees) {
    return degrees * M_PI / 180.0;
}

__device__ __forceinline__ double haversine(double lon1, double lat1, double lon2, double lat2) {
    // 将经纬度从度转换为弧度
    lon1 = degreesToRadians(lon1);
    lat1 = degreesToRadians(lat1);
    lon2 = degreesToRadians(lon2);
    lat2 = degreesToRadians(lat2);

    // 差值
    double dlon = lon2 - lon1;
    double dlat = lat2 - lat1;

    // 哈弗赛因公式
    double a = std::sin(dlat / 2) * std::sin(dlat / 2) +
               std::cos(lat1) * std::cos(lat2) *
               std::sin(dlon / 2) * std::sin(dlon / 2);
    double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));

    // 计算距离
    double distance = EARTH_RADIUS_KM * c;

    return distance;
}


__device__ __forceinline__ double gpu_max_distance(Point &p, box &bx){
	Point q;
	q.x = (abs(p.x-bx.low[0]) < abs(p.x-bx.high[0])) ? bx.high[0] : bx.low[0];
	q.y = (abs(p.y-bx.low[1]) < abs(p.y-bx.high[1])) ? bx.high[1] : bx.low[1];

	// printf("%lf %lf %lf %lf %lf\n", p.x, p.y, q.x, q.y, haversine(p.x, p.y, q.x, q.y));

	return haversine(p.x, p.y, q.x, q.y);
}

__device__ __forceinline__ double gpu_point_to_segment_distance(const Point &p, const Point &p1, const Point &p2){
    double A = p.x - p1.x;
    double B = p.y - p1.y;
    double C = p2.x - p1.x;
    double D = p2.y - p1.y;

    double dot = A * C + B * D;
    double len_sq = C * C + D * D;
    double param = -1;
    if (len_sq != 0) // in case of 0 length line
        param = dot / len_sq;

    double xx, yy;

    if (param < 0) {
        xx = p1.x;
        yy = p1.y;
    } else if (param > 1) {
        xx = p2.x;
        yy = p2.y;
    } else {
        xx = p1.x + param * C;
        yy = p1.y + param * D;
    }

	return haversine(p.x, p.y, xx, yy);
}