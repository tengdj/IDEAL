/*
 * partition.cpp
 *
 *  Created on: May 28, 2022
 *      Author: teng
 *
 *
 *  this tool helps to partition the space into N*M tiles based on the distribution of the points
 *
 */

#include "../index/RTree.h"
#include <queue>
#include <fstream>
#include "../geometry/MyPolygon.h"
#include "../geometry/query_context.h"
#include <vector>
#include "../util/util.h"

using namespace std;

int dimx = 20;
int dimy = 20;

bool comparePointX(Point p1, Point p2)
{
    return (p1.x < p2.x);
}

bool comparePointY(Point p1, Point p2)
{
    return (p1.y < p2.y);
}

int main(int argc, char** argv) {
	query_context global_ctx;
	global_ctx = get_parameters(argc, argv);
	global_ctx.load_points();

	int num = global_ctx.target_num;
	//int num = 10000000;
	dimx = 100;
	dimy = 100;
	vector<Point> points((Point *)global_ctx.points, ((Point *)global_ctx.points)+num);

	struct timeval start = get_cur_time();
	sort(points.begin(), points.end(), comparePointX);
	logt("sort %d points", start, num);

	for(int x=0;x<dimx;x++){
		int begin = (num/dimx)*x;
		int end = (x+1)*(num/dimx);
		if(x==dimx-1){
			end = num;
		}
		sort(points.begin()+begin, points.begin()+end, comparePointY);
		logt("sort %d to %d (%d)", start, begin, end, end-begin);
	}

	Pixel parts[dimx][dimy];
	for(int x=0;x<dimx;x++){
		int size = num/dimx;
		if(x==dimx-1){
			size = num - x*(num/dimx);
		}
		log("%d",size);
		for(int y=0;y<dimy;y++){
			int begin = x*(num/dimx)+y*(size/dimy);
			int end = x*(num/dimx)+(y+1)*(size/dimy);
			end = min(end,num);
			for(int t = begin;t<end;t++){
				//points[pids[t]].print();
				parts[x][y].update(points[t]);
			}

			int pid = x*dimy+y;
			char path[256];
			sprintf(path, "part_points/%d.dat",pid);
			ofstream os;
			os.open(path, ios::out | ios::binary |ios::trunc);
			assert(os.is_open());
			Point *pa = &points[0];
			os.write((char *)(pa+begin), sizeof(Point)*(end-begin));
			os.close();
		}
	}

	vector<Pixel *> pixels;
	for(int x=0;x<dimx;x++){
		for(int y=0;y<dimy;y++){
			pixels.push_back(&parts[x][y]);
		}
	}
	print_boxes(pixels);

	return 0;
}



