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

#define VECTORS_PER_WORKITEM 16

class OpenCLKernel: public AlignmentKernel {

public:

	OpenCLKernel() {
		bool exception = false;

		Parameters.has_key("score_match") ?
				scoreMatch = Parameters.param_int("score_match") : exception =
						true;
		Parameters.has_key("score_mismatch") ?
				scoreMismatch = Parameters.param_int("score_mismatch") :
				exception = true;
		Parameters.has_key("score_gap_read") ?
				scoreGapRead = Parameters.param_int("score_gap_read") :
				exception = true;
		Parameters.has_key("score_gap_ref") ?
				scoreGapRef = Parameters.param_int("score_gap_ref") :
				exception = true;
		Parameters.has_key("read_length") ?
				readLength = Parameters.param_int("read_length") : exception =
						true;
		Parameters.has_key("ref_length") ?
				refLength = Parameters.param_int("ref_length") : exception =
						true;

		if (exception) {
			throw "Cannot instantiate Kernel. Lacking parameters";
		}

		alnLength = refLength + readLength;
		matrix_size = (readLength + 1) * (refLength + 1);

		initialize_opencl_environment();
	}
	~OpenCLKernel() {
		if (host_reads != 0 && host_refs != 0 && host_results != 0
				&& host_indices != 0 && host_matrix != 0) {
			delete[] host_reads;
			delete[] host_refs;
			delete[] host_results;
			delete[] host_indices;
			delete[] host_matrix;
		}
	}
	void score_alignments(int const & opt, int const & aln_number,
			char const * const * const reads, char const * const * const refs,
			short * const scores);
	void compute_alignments(int const & opt, int const & aln_number,
			char const * const * const reads, char const * const * const refs,
			Alignment * const alignments);

private:

	cl::Device device;
	cl::Context context;
	cl::Program program;
	cl::CommandQueue queue;
	cl::Kernel kernel;

	cl::Device setup_opencl_device(cl_device_type const & device_type);
	std::vector<cl::Device> fission_opencl_device(cl::Device const & device);
	cl::Context setup_context(cl::Device const & device);
	cl::Program setup_program(cl::Context const & context);
	cl::CommandQueue setup_queue(cl::Context const & context,
			cl::Device const & device);
	cl::Kernel setup_kernel(cl::Program const & program,
			char const * kernel_name);
	size_t calculate_batch_size_from_memory(cl::Kernel const & kernel,
			cl::Device const & device, int const & opt);
	void partition_load(int const & aln_number, size_t const & batch_size,
			size_t & batch_num, size_t & overhang);

	void initialize_opencl_environment();
	void init_host_memory(size_t const & batch_size, bool const & score);

	int readLength;
	int refLength;
	int alnLength;
	int matrix_size;

	short scoreMatch;
	short scoreMismatch;
	short scoreGapRead;
	short scoreGapRef;

	char * host_reads = 0;
	char * host_refs = 0;
	char * host_results = 0;
	short * host_indices = 0;
	short * host_matrix = 0;

	char const * const score_smith_waterman_kernel = "score_alignment_smith_waterman";
	char const * const align_smith_waterman_kernel = "calc_alignment_smith_waterman";

	char const * const score_needleman_wunsch_kernel = "score_alignment_needleman_wunsch";
	char const * const align_needleman_wunsch_kernel = "calc_alignment_needleman_wunsch";

	char const * cl_error_to_string(cl_int ci_error_num);

	void check_opencl_success(char const * msg, cl_int ci_error_num);
};

#endif /* OPENCLKERNEL_H */
