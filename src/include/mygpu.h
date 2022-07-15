/*
 * mygpu.h
 *
 *  Created on: Dec 9, 2019
 *      Author: teng
 */

#ifndef MYGPU_H_
#define MYGPU_H_

//#ifdef USE_GPU

#include <vector>
#include <pthread.h>
#include <map>

#include "util.h"

using namespace std;

class gpu_info{
	pthread_mutex_t lk;
	map<size_t, size_t> freed;
public:
	int device_id;
	size_t mem_size;
	bool busy;

	double *degree_per_kilometer;
	char *d_data = NULL;
	size_t data_size = 0;
	size_t allocated_begin = 0;
	size_t allocated_end = 0;

	void lock();
	void unlock();
	gpu_info(int id, size_t mem_size);
	~gpu_info();
	char *allocate(size_t ss);
	void free(char *d_addr, size_t offset);
};

extern vector<gpu_info *> gpus;

void init_gpus();
void release_gpus();
vector<gpu_info *> get_gpus();
void print_gpus();

//#endif


#endif /* MYGPU_H_ */
