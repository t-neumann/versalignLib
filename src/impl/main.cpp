/*
 * main.cpp
 *
 *  Created on: May 20, 2016
 *      Author: tobias.neumann
 */

#include "AlignmentKernel.h"
#include "CustomParameters.h"
#include "CustomLogger.h"

#include "util/versalignUtil.h"

#include <iostream>
#include <fstream>
#include <map>

#define READ_PAD '\0'
#define REF_PAD '\0'

using std::cout;
using std::string;
using std::endl;

//--------------//
// Defines
//--------------//

#ifndef NDEBUG
const string libPathDir = "../bin/Versalign-0.1.0-debug/lib/";
#else
const string libPathDir = "../bin/Versalign-0.1.0/lib/";
#endif

#ifdef __APPLE__
const string libSuffix(".dylib");
#else
const string libSuffix (".so");
#endif

//--------------//
// Helpers
//--------------//

AlignmentKernel * get_kernel(int const & dll);
void clear_kernel(AlignmentKernel * kernel, char const * kernel_name);

void time_kernel(char const * kernel_name, char const * * const reads,
		char const * * const refs, size_t const & seqNumber,
		bool const & score = true);

inline bool lib_exists(const std::string& name) {
	std::ifstream f(name.c_str());
	return f.good();
}

//--------------//
// Globals
//--------------//

std::map<std::string, std::string> kernel_map =
		{ { "default", libPathDir + "libDefaultKernel" + libSuffix }, { "SSE",
				libPathDir + "libSSEKernel" + libSuffix }, { "OpenCL",
				libPathDir + "libOpenCLKernel" + libSuffix } };

int const num_loops = 100;

int const trials = 7;
int const threads[trials] = { 1, 2, 4, 8, 16, 32, 64 };

CustomParameters parameters;
CustomLogger logger;

int main(int argc, char *argv[]) {

# ifdef _OPENMP
	logger.log(0, "MAIN", "Compiled by an OpenMP-compliant implementation.");
# endif

	if (check_avx2_support()) {
		logger.log(0, "MAIN", "AVX2 boost detected.");
		kernel_map["AVX"] = libPathDir + "libAVXKernel" + libSuffix;
	}

	FastaProvider read_provider;
	std::vector<const char *> readsVec = read_provider.parse_fasta(
			"../testset/reads.fa");
	char const * * const reads = &readsVec[0];
	std::vector<const char *> refsVec = read_provider.parse_fasta(
			"../testset/refs.fa");
	char const * * const refs = &refsVec[0];

	if (readsVec.size() != refsVec.size()) {

		logger.log(3, "MAIN",
				std::string(
						"Unequal sizes of reads and ref set ("
								+ std::to_string(readsVec.size()) + " vs "
								+ std::to_string(refsVec.size()) + ").").c_str());

		return -1;
	}

	size_t const seqNumber = readsVec.size();

	size_t const max_read_length = pad(reads, seqNumber, READ_PAD);
	size_t const max_ref_length = pad(refs, seqNumber, REF_PAD);

	parameters.read_length = max_read_length;
	parameters.ref_length = max_ref_length;
	parameters.num_threads = 10;

	//////////////////////////////
	// Example calls
	//////////////////////////////

	char const * kernel_name = "OpenCL";

	int dll = DLL_init(kernel_map[kernel_name].c_str(), &parameters, &logger);
	AlignmentKernel * kernel = get_kernel(dll);

	short * scores = new short[seqNumber]();
	Alignment * alignments = new Alignment[seqNumber]();

	//////////////////////////////
	// Smith Waterman-Mode
	//////////////////////////////

	int mode = 0;

	kernel->score_alignments(mode, seqNumber, reads, refs, scores);

	std::ofstream outfile;
	outfile.open("scores_smith_waterman.txt");

	for (int i = 0; i < seqNumber; ++i) {
		outfile << reads[i] << "\t" << scores[i] << std::endl;

	}

	outfile.close();

	kernel->compute_alignments(mode, seqNumber, reads, refs, alignments);

	outfile.open("alignments_smith_waterman.txt");

	for (int i = 0; i < seqNumber; ++i) {
		outfile << alignments[i].read + alignments[i].readStart << std::endl;
		outfile << alignments[i].ref + alignments[i].refStart << std::endl
				<< std::endl;
	}

	outfile.close();

	clear_kernel(kernel, kernel_name);

	//////////////////////////////
	// Needelman Wunsch-Mode
	//////////////////////////////

	mode = 1;

	kernel = get_kernel(dll);

	kernel->score_alignments(mode, seqNumber, reads, refs, scores);

	outfile.open("scores_needleman_wunsch.txt");

	for (int i = 0; i < seqNumber; ++i) {
		outfile << reads[i] << "\t" << scores[i] << std::endl;

	}

	outfile.close();

	kernel->compute_alignments(mode, seqNumber, reads, refs, alignments);

	outfile.open("alignments_needleman_wunsch.txt");

	for (int i = 0; i < seqNumber; ++i) {
		outfile << alignments[i].read + alignments[i].readStart << std::endl;
		outfile << alignments[i].ref + alignments[i].refStart << std::endl
				<< std::endl;
	}

	outfile.close();

	clear_kernel(kernel, kernel_name);

	//////////////////////////////
	// Timing of different Kernels
	//////////////////////////////

	std::cout << "Threads";

	for (int trial = 0; trial < trials; ++trial) {

		std::cout << "\t" << threads[trial];

	}

	std::cout << std::endl;

	for (std::map<std::string, std::string>::iterator i = kernel_map.begin();
			i != kernel_map.end(); i++) {

		time_kernel(i->first.c_str(), reads, refs, seqNumber, false);

	}

	return 0;
}

void clear_kernel(AlignmentKernel * kernel, char const * kernel_name) {
	int dll = DLL_init(kernel_map[kernel_name].c_str(), &parameters, &logger);

	fp_delete_alignment_kernel delete_alignment_kernel =
			(fp_delete_alignment_kernel) DLL_function_retreival(dll,
					"delete_alignment_kernel");

	delete_alignment_kernel(kernel);
}

AlignmentKernel * get_kernel(int const & dll) {

	fp_load_alignment_kernel load_alignment_kernel =
			(fp_load_alignment_kernel) DLL_function_retreival(dll,
					"spawn_alignment_kernel");

	AlignmentKernel * kernel = 0;

	kernel = load_alignment_kernel();

	return kernel;
}

void time_kernel(char const * kernel_name, char const * * const reads,
		char const * * const refs, size_t const & seqNumber,
		bool const & score) {

	int dll = DLL_init(kernel_map[kernel_name].c_str(), &parameters, &logger);
	AlignmentKernel * kernel = get_kernel(dll);

	fp_set_parameters set_parameters =
			(fp_set_parameters) DLL_function_retreival(dll, "set_parameters");

	short * scores = new short[seqNumber]();
	Alignment * alignments = new Alignment[seqNumber]();

	std::cout << kernel_name;

	Timer timer;

	for (int trial = 0; trial < trials; ++trial) {

		int const num_threads = threads[trial];

		parameters.num_threads = num_threads;

		set_parameters(&parameters);

		AlignmentKernel * kernel = get_kernel(dll);

		if (score) {
			timer.start();

			for (int rep = 0; rep < num_loops; ++rep) {

				kernel->score_alignments(0, seqNumber, reads, refs, scores);

			}

			timer.stop();
		} else {
			timer.start();

			for (int rep = 0; rep < num_loops; ++rep) {

				kernel->compute_alignments(0, seqNumber, reads, refs,
						alignments);

			}

			timer.stop();
		}

		clear_kernel(kernel, kernel_name);

		std::cout << "\t" << timer.getElapsedTimeInMicroSec() / num_loops;
	}
	std::cout << std::endl;
}
