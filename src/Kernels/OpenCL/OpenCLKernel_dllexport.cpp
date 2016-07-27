/*
 * OpenCLKernel_dllexport.cpp
 *
 *  Created on: Jul 25, 2016
 *      Author: elpresidente
 */

#include "OpenCLKernel.h"

#ifdef _WIN32
#define dllexport  __declspec(dllexport)
#else
#define dllexport
#endif

extern "C" dllexport AlignmentKernel * spawn_alignment_kernel() {

	AlignmentKernel * kernel = new OpenCLKernel();

	return kernel;
}

AlignmentParameters * _parameters = 0;

extern "C" dllexport void set_parameters(AlignmentParameters * parameters) {
	_parameters = parameters;
}

AlignmentLogger * _logger = 0;

extern "C" dllexport void set_logger(AlignmentLogger * logger) {
	_logger = logger;
}

extern "C" dllexport void delete_alignment_kernel(OpenCLKernel * instance) {
	if (instance != 0) {
		delete instance;
		instance = 0;
	}
}
