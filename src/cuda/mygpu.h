/*
 * mygpu.h
 *
 *  Created on: Dec 9, 2019
 *      Author: teng
 */

#ifndef MYGPU_H_
#define MYGPU_H_

#include <pthread.h>
#include <vector>
using namespace std;


class gpu_info{
public:
	int device_id;
	size_t mem_size;
	bool busy;
	pthread_mutex_t lock;
	double *d_data = NULL;
	size_t data_size = 0;
	double *source_data = NULL;
	size_t source_size = 0;
	int *result = NULL;
	size_t result_size = 0;
	uint *offset_size = NULL;
	size_t os_size = 0;


	void init();
	~gpu_info();
	double *get_source(size_t ss);
	double *get_data(size_t ds);
	int *get_result(size_t rs);
	uint *get_os(size_t os);

};

vector<gpu_info *> get_gpus();
void print_gpus();


#endif /* MYGPU_H_ */
