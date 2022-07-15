/*
 * cuda_util.cuh
 *
 *  Created on: Jul 14, 2022
 *      Author: teng
 */

#ifndef SRC_CUDA_CUDA_UTIL_CUH_
#define SRC_CUDA_CUDA_UTIL_CUH_

#include <cuda.h>
#include "mygpu.h"



void check_execution();

#define CUDA_SAFE_CALL(call) 										  	  \
	do {																  \
		cudaError_t err = call;											  \
		if (cudaSuccess != err) {										  \
			fprintf (stderr, "Cuda error in file '%s' in line %i : %s.\n",\
					__FILE__, __LINE__, cudaGetErrorString(err) );	      \
			exit(EXIT_FAILURE);											  \
		}																  \
	} while (0);

inline void check_execution(){
	cudaError_t err = cudaGetLastError();
	if (err != cudaSuccess){
		log(cudaGetErrorString(err));
	}
}

#endif /* SRC_CUDA_CUDA_UTIL_CUH_ */
