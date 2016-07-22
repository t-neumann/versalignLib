/*
 * OpenCLKernel.h
 *
 *  Created on: Jul 22, 2016
 *      Author: tobias.neumann
 */

#ifndef OPENCLKERNEL_H
#define OPENCLKERNEL_H

#include "AlignmentKernel.h"
#include "AlignmentParameters.h"

#include <CL/cl.hpp>

class OpenCLKernel : public AlignmentKernel {

	public:

		OpenCLKernel() {
			bool exception = false;

			Parameters.has_key("score_match") ? scoreMatch = Parameters.param_int("score_match") : exception = true;
			Parameters.has_key("score_mismatch") ? scoreMismatch = Parameters.param_int("score_mismatch") : exception = true;
			Parameters.has_key("score_gap_read") ? scoreGapRead = Parameters.param_int("score_gap_read") : exception = true;
			Parameters.has_key("score_gap_ref") ? scoreGapRef = Parameters.param_int("score_gap_ref") : exception = true;
			Parameters.has_key("read_length") ? readLength = Parameters.param_int("read_length") : exception = true;
			Parameters.has_key("ref_length") ? refLength = Parameters.param_int("ref_length") : exception = true;

			alnLength = refLength + readLength;

			if (exception) {
				throw "Cannot instantiate Kernel. Lacking parameters";
			}

			device = setup_opencl_device(CL_DEVICE_TYPE_CPU);
			context = setup_context(device);
		}
		~OpenCLKernel() {}
		void score_alignments(int const & opt, int const & aln_number, char const * const * const reads,
				char const * const * const refs, short * const scores);
		void compute_alignments(int const & opt, int const & aln_number, char const * const * const reads,
				char const * const * const refs, Alignment * const alignments);

	private:

		cl::Device device;
		cl::Context context;
		cl::Program program;

		cl::Device setup_opencl_device(cl_device_type const & device_type);
		cl::Context setup_context(cl::Device const & device);
		cl::Program setup_program(cl::Context const & context);
		cl::CommandQueue setup_queue(cl::Context const & context, cl::Device const & device);

		int readLength;
		int refLength;
		int alnLength;

		short scoreMatch;
		short scoreMismatch;
		short scoreGapRead;
		short scoreGapRef;

		char const * cl_error_to_string(cl_int ci_error_num);

		void check_opencl_success(char const * msg, cl_int ci_error_num);
};



#endif /* OPENCLKERNEL_H */
