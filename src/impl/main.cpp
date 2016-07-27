/*
 * main.cpp
 *
 *  Created on: May 20, 2016
 *      Author: tobias.neumann
 */

#include "AlignmentKernel.h"
#include "CustomParameters.h"
#include "CustomLogger.h"
#include "dll_handling.h"
//#include "Kernels/plain/SWKernel.h"
//#include "Kernels/AVX-SSE/SSEKernel.h"
#include "Kernels/AVX-SSE/AVXKernel.h"
#include "util/versalignUtil.h"
//#include "Kernels/OpenCL/ocl_testing.h"
#include "Timer/Timer.h"
#include "AVXSupportChecker.h"
#include "FastaProvider.h"
//#include "Kernels/OpenCL/OpenCLKernel.h"
#include <iostream>
#include <fstream>

#define READ_PAD '\0'
#define REF_PAD '\0'

using std::cout;
using std::string;
using std::endl;

CustomParameters parameters;
CustomLogger logger;

inline bool lib_exists(const std::string& name) {
	std::ifstream f(name.c_str());
	return f.good();
}

int main(int argc, char *argv[]) {

//		char const
//		* reads[] =
//		{
//				"AGGGGGGA",
//				"AATTTTGCC",
//				"TTTTTAA",
//				"ATAGATAGATAG",
//				"AGCAGTAC",
//				"AGAGAGAG",
//				"",
//				"ATATATAT",
//				"AGG",
//				"TAA",
//				"AGAG",
//				"AGCAGTAG",
//				"TATAC",
//				"AGGAAGAG",
//				"TT",
//				"CCCCC",
//				"AGGGGGGA",
//				"AGGGGGGA",
//				"AATTTTGCC",
//				"TTTTTAA",
//				"ATAGATAGATAG",
//				"AGCAGTAC",
//				"AGAGAGAG",
//				"",
//				"ATATATAT",
//				"AGG",
//				"TAA",
//				"AGAG",
//				"AGCAGTAG",
//				"TATAC",
//				"AGGAAGAG",
//				"TT",
//				"CCCCC",
//				"AGGGGGGA",
//				"AGGGGGGA",
//				"AATTTTGCC",
//				"TTTTTAA",
//				"ATAGATAGATAG",
//				"AGCAGTAC",
//				"AGAGAGAG",
//				"",
//				"ATATATAT",
//				"AGG",
//				"TAA",
//				"AGAG",
//				"AGCAGTAG",
//				"TATAC",
//				"AGGAAGAG",
//				"TT",
//				"CCCCC",
//				"AGGGGGGA",
//				"AGGGGGGA",
//				"AATTTTGCC",
//				"TTTTTAA",
//				"ATAGATAGATAG",
//				"AGCAGTAC",
//				"AGAGAGAG",
//				"",
//				"ATATATAT",
//				"AGG",
//				"TAA",
//				"AGAG",
//				"AGCAGTAG",
//				"TATAC",
//				"AGGAAGAG",
//				"TT",
//				"CCCCC",
//				"AGGGGGGA",
//				"AGGGGGGA",
//				"AATTTTGCC",
//				"TTTTTAA",
//				"ATAGATAGATAG",
//				"AGCAGTAC",
//				"AGAGAGAG",
//				"",
//				"ATATATAT",
//				"AGG",
//				"TAA",
//				"AGAG",
//				"AGCAGTAG",
//				"TATAC",
//				"AGGAAGAG",
//				"TT",
//				"CCCCC",
//				"AGGGGGGA"
//		};
//
//		char const
//		* refs[] =
//		{
//				"AGCAGATGAGGGCGGATAGC",
//				"TTTTGCCAACGCATGGCAGA",
//				"ATGACGACGCAGTGCTTTTT",
//				"GCAATAGAATAGATAGTGGT",
//				"TAGCATCAAGCAGTACTACA",
//				"TCTCTCTCTCTCTCTCTCTC",
//				"AGCAGATGACATGCATGCAA",
//				"",
//				"AGCAGATGAGGGCGGATAGC",
//				"TTTTGCCAACGCATGGCAGA",
//				"ATGACGACGCAGTGCTTTTT",
//				"GCAATAGAATAGATAGTGGT",
//				"TAGCATCAAGCAGTACTACA",
//				"TCTCTCTCTCTCTCTCTCTC",
//				"AGCAGATGACATGCATGCAA",
//				"",
//				"AGCAGATGAGGGCGGATAGC",
//				"AGCAGATGAGGGCGGATAGC",
//				"TTTTGCCAACGCATGGCAGA",
//				"ATGACGACGCAGTGCTTTTT",
//				"GCAATAGAATAGATAGTGGT",
//				"TAGCATCAAGCAGTACTACA",
//				"TCTCTCTCTCTCTCTCTCTC",
//				"AGCAGATGACATGCATGCAA",
//				"",
//				"AGCAGATGAGGGCGGATAGC",
//				"TTTTGCCAACGCATGGCAGA",
//				"ATGACGACGCAGTGCTTTTT",
//				"GCAATAGAATAGATAGTGGT",
//				"TAGCATCAAGCAGTACTACA",
//				"TCTCTCTCTCTCTCTCTCTC",
//				"AGCAGATGACATGCATGCAA",
//				"",
//				"AGCAGATGAGGGCGGATAGC",
//				"AGCAGATGAGGGCGGATAGC",
//				"TTTTGCCAACGCATGGCAGA",
//				"ATGACGACGCAGTGCTTTTT",
//				"GCAATAGAATAGATAGTGGT",
//				"TAGCATCAAGCAGTACTACA",
//				"TCTCTCTCTCTCTCTCTCTC",
//				"AGCAGATGACATGCATGCAA",
//				"",
//				"AGCAGATGAGGGCGGATAGC",
//				"TTTTGCCAACGCATGGCAGA",
//				"ATGACGACGCAGTGCTTTTT",
//				"GCAATAGAATAGATAGTGGT",
//				"TAGCATCAAGCAGTACTACA",
//				"TCTCTCTCTCTCTCTCTCTC",
//				"AGCAGATGACATGCATGCAA",
//				"",
//				"AGCAGATGAGGGCGGATAGC",
//				"AGCAGATGAGGGCGGATAGC",
//				"TTTTGCCAACGCATGGCAGA",
//				"ATGACGACGCAGTGCTTTTT",
//				"GCAATAGAATAGATAGTGGT",
//				"TAGCATCAAGCAGTACTACA",
//				"TCTCTCTCTCTCTCTCTCTC",
//				"AGCAGATGACATGCATGCAA",
//				"",
//				"AGCAGATGAGGGCGGATAGC",
//				"TTTTGCCAACGCATGGCAGA",
//				"ATGACGACGCAGTGCTTTTT",
//				"GCAATAGAATAGATAGTGGT",
//				"TAGCATCAAGCAGTACTACA",
//				"TCTCTCTCTCTCTCTCTCTC",
//				"AGCAGATGACATGCATGCAA",
//				"",
//				"AGCAGATGAGGGCGGATAGC",
//				"AGCAGATGAGGGCGGATAGC",
//				"TTTTGCCAACGCATGGCAGA",
//				"ATGACGACGCAGTGCTTTTT",
//				"GCAATAGAATAGATAGTGGT",
//				"TAGCATCAAGCAGTACTACA",
//				"TCTCTCTCTCTCTCTCTCTC",
//				"AGCAGATGACATGCATGCAA",
//				"",
//				"AGCAGATGAGGGCGGATAGC",
//				"TTTTGCCAACGCATGGCAGA",
//				"ATGACGACGCAGTGCTTTTT",
//				"GCAATAGAATAGATAGTGGT",
//				"TAGCATCAAGCAGTACTACA",
//				"TCTCTCTCTCTCTCTCTCTC",
//				"AGCAGATGACATGCATGCAA",
//				"",
//				"AGCAGATGAGGGCGGATAGC"
//		};
//
//		int seqNumber = 17;

	int const N = 100;

	FastaProvider read_provider;
	std::vector<const char *> readsVec = read_provider.parse_fasta(
			"../testset/reads.fa");
	char const * * const reads = &readsVec[0];
	std::vector<const char *> refsVec = read_provider.parse_fasta(
			"../testset/refs.fa");
	char const * * const refs = &refsVec[0];

	if (readsVec.size() != refsVec.size()) {

		logger.log(3, "MAIN", std::string("Unequal sizes of reads and ref set (" + std::to_string(readsVec.size()) + " vs " + std::to_string(refsVec.size()) + ").").c_str());

		return -1;
	}

	int seqNumber = readsVec.size();

	size_t max_read_length = pad(reads, seqNumber, READ_PAD);

	size_t max_ref_length = pad(refs, seqNumber, REF_PAD);

	parameters.read_length = max_read_length;
	parameters.ref_length = max_ref_length;
	parameters.num_threads = 10;

#ifndef NDEBUG
	const string libPathDir = "../bin/Versalign-0.1.0-debug/lib/";
#else
	const string libPathDir = "../bin/Versalign-0.1.0/lib/";
#endif

#ifdef __APPLE__
	const string libSuffix("dylib");
#else
	const string libSuffix ("so");
#endif

# ifdef _OPENMP
	logger.log(0, "MAIN", "Compiled by an OpenMP-compliant implementation.");

# endif

	string libPath(libPathDir + "libDefaultKernel." + libSuffix);

	if (check_avx2_support() && lib_exists(libPath)) {
		logger.log(0, "MAIN", "AVX2 boost detected.");
//		libPath = libPathDir + "libAVXKernel." + libSuffix;
	} else {
		logger.log(0, "MAIN", "SSE4 mode activated.");
//		libPath = libPathDir + "libSSEKernel." + libSuffix;
	}

//	libPath = libPathDir + "libOpenCLKernel." + libSuffix;

	int dll = DLL_init(libPath.c_str(), &parameters, &logger);

	//OpenCLKernel * ocl = new OpenCLKernel();

	fp_load_alignment_kernel load_alignment_kernel =
			(fp_load_alignment_kernel) DLL_function_retreival(dll,
					"spawn_alignment_kernel");
	fp_delete_alignment_kernel delete_alignment_kernel =
			(fp_delete_alignment_kernel) DLL_function_retreival(dll,
					"delete_alignment_kernel");
	fp_set_parameters set_parameters =
			(fp_set_parameters) DLL_function_retreival(dll, "set_parameters");

	AlignmentKernel * plain_kernel = 0;

	plain_kernel = load_alignment_kernel();

	logger.log(0, "MAIN", "Alignment kernel instantiated.");
	logger.log(0, "MAIN", std::string("Running with " + std::to_string(parameters.num_threads) + " threads.").c_str());

	short * scores = new short[seqNumber]();
	Alignment * alignments = new Alignment[seqNumber]();

//	run_ocl_test(reads, refs, seqNumber, max_read_length, max_ref_length);

//	ocl->score_alignments(0, seqNumber, reads, refs, scores);
//	ocl->compute_alignments(1, seqNumber, reads, refs, alignments);

	plain_kernel->compute_alignments(1, seqNumber, reads, refs, alignments);
//	plain_kernel->score_alignments(0, seqNumber, reads, refs, scores);

//	for (int i = 0; i < seqNumber; ++i) {
//		std::cout << reads[i] << std::endl;
//		std::cout << refs[i] << std::endl;
//		std::cout << scores[i] << std::endl;
//	}

//	for (int i = 0; i < seqNumber; ++i) {
//		std::cout << "Read: " << reads[i] << std::endl;
//		std::cout << "Ref: " << refs[i] << std::endl;
//		std::cout << "==================" << std::endl << "\"";
//		std::cout << alignments[i].read + alignments[i].readStart;
//		std::cout << "\"" << std::endl << "\"";
//		std::cout << alignments[i].ref + alignments[i].refStart;
//		std::cout << "\"" << std::endl << "==================" << std::endl;
//		std::cout << "Score: " << scores[i] << std::endl << std::endl;
//	}

	//plain_kernel->score_alignments(0, seqNumber, reads,refs, scores);

	return 0;

	Timer timer;

	timer.start();

	for (int rep = 0; rep < N; ++rep) {

		//plain_kernel->compute_alignments(0, seqNumber, reads, refs, alignments);
		plain_kernel->score_alignments(0, seqNumber, reads, refs, scores);

	}

	timer.stop();

	cout << "Alignment took " << timer.getElapsedTimeInMicroSec() / N << " ms"
			<< endl;

	parameters.num_threads = 4;
	set_parameters(&parameters);

	std::cout << "Alignment kernel instantiated.\n";
	std::cout << "Running with " << parameters.num_threads << " threads.\n";

	timer.start();

	for (int rep = 0; rep < N; ++rep) {

		//plain_kernel->compute_alignments(0, seqNumber, reads, refs, alignments);
		plain_kernel->score_alignments(0, seqNumber, reads, refs, scores);

	}

	timer.stop();
	cout << "Alignment took " << timer.getElapsedTimeInMicroSec() / N << " ms"
			<< endl;

	parameters.num_threads = 8;
	set_parameters(&parameters);

	std::cout << "Alignment kernel instantiated.\n";
	std::cout << "Running with " << parameters.num_threads << " threads.\n";

	timer.start();

	for (int rep = 0; rep < N; ++rep) {

		//plain_kernel->compute_alignments(0, seqNumber, reads, refs, alignments);
		plain_kernel->score_alignments(0, seqNumber, reads, refs, scores);

	}

	timer.stop();
	cout << "Alignment took " << timer.getElapsedTimeInMicroSec() / N << " ms"
			<< endl;

//	for (int i = 0; i < seqNumber; ++i) {
//		std::cout << "Read: " << reads[i] << std::endl;
//		std::cout << "Ref: " << refs[i] << std::endl;
//		std::cout << "==================" << std::endl << "\"";
//		std::cout << alignments[i].read + alignments[i].readStart;
//		std::cout << "\"" << std::endl << "\"";
//		std::cout << alignments[i].ref + alignments[i].refStart;
//		std::cout << "\"" << std::endl << "==================" << std::endl;
//		std::cout << "Score: " << scores[i] << std::endl << std::endl;
//	}

	delete_alignment_kernel(plain_kernel);

	return 0;
}

