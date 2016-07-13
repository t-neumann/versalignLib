/*
 * SWKernel_dllexport.cpp
 *
 *  Created on: Jul 13, 2016
 *      Author: Tobias Neumann
 *       Email: tobias.neumann.at@gmail.com
 */

#include "SWKernel.h"
#include "AlignmentKernel.h"

#ifdef _WIN32
#define dllexport  __declspec(dllexport)
#else
#define dllexport
#endif

extern "C" dllexport AlignmentKernel * spawn_alignment_kernel() {

	AlignmentKernel * kernel = new SWKernel();

	return kernel;
}


extern "C" dllexport void delete_alignment_kernel(SWKernel * instance) {
	if (instance != 0) {
			delete instance;
			instance = 0;
		}
}
