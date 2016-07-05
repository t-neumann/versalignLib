/*
 * ocl_testing.h
 *
 *  Created on: Jun 29, 2016
 *      Author: tobias.neumann
 */

#ifndef SRC_KERNELS_OPENCL_OCL_TESTING_H_
#define SRC_KERNELS_OPENCL_OCL_TESTING_H_

#include <CL/cl.hpp>
#include <string>
#include <vector>
#include <iostream>

void run_ocl_test(char const * const * const reads, char const * const * const refs, int const & seqNumber, size_t const & max_read_length, size_t const & max_ref_length);


#endif /* SRC_KERNELS_OPENCL_OCL_TESTING_H_ */
