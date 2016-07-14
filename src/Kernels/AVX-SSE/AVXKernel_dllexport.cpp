/*
 * AVXKernel_dllexport.cpp
 *
 *  Created on: Jul 14, 2016
 *      Author: tobias.neumann
 */

#include "AVXKernel.h"
#include "AlignmentKernel.h"

#ifdef _WIN32
#define dllexport  __declspec(dllexport)
#else
#define dllexport
#endif

extern "C" dllexport AlignmentKernel * spawn_alignment_kernel() {

	AlignmentKernel * kernel = new AVXKernel();

	return kernel;
}


AlignmentParameters * _parameters = 0;

extern "C" dllexport void set_parameters(AlignmentParameters * parameters) {
	_parameters = parameters;
}


extern "C" dllexport void delete_alignment_kernel(AVXKernel * instance) {
	if (instance != 0) {
		delete instance;
		instance = 0;
	}
}


