/*
 * main.cpp
 *
 *  Created on: May 20, 2016
 *      Author: tobias.neumann
 */

#include "Kernels/plain/SWKernel.h"
#include "Kernels/SSE/SSEKernel.h"
#include "util/versalignUtil.h"
#include "Kernels/OpenCL/ocl_testing.h"
#include "Timer/Timer.h"
#include <iostream>

#define READ_PAD '\0'
#define REF_PAD '\0'

using std::cout;
using std::string;
using std::endl;

//#ifdef __APPLE__
//    #include "OpenCL/opencl.h"
//#else
//    #include <CL/cl.h>
//#endif

int main(int argc, char *argv[]) {

	run_ocl_test();

	return 0;

	/*char const
	 * reads[] =
						{
								"CACACCCACACACCACACCACACACCAGACCCACACCCACACACACACATCCTAAGACTGCCCTAAAACAGCCCTAATCTAACCCTGGCCAACCTGTCTCT",
								"CACACCCACACACCACACCACACACCAGACCCACACCCACACACACACATCCTAAGACTGCCCTAAAACAGCCCTAATCTAACCCTGGCCAACCTGTCTCT",
								"CACACCCACACACCACACCACACACCAGACCCACACCCACACACACACATCCTAAGACTGCCCTAAAACAGCCCTAATCTAACCCTGGCCAACCTGTCTCT",
								"CACACCCACACACCACACCACACACCAGACCCACACCCACACACACACATCCTAAGACTGCCCTAAAACAGCCCTAATCTAACCCTGGCCAACCTGTCTCT",
								"CACACCCACACACCACACCACACACCAGACCCACACCCACACACACACATCCTAAGACTGCCCTAAAACAGCCCTAATCTAACCCTGGCCAACCTGTCTCT",
								"GGGTAAGTTGAGAGACAGGTTGGACAGGGTTAGATTAGGGCTGTGTTAGGGTAGTGTTAGGATGTGTGTGTGTGGGTGTGGTGTGGTGTGTGGTGTGGTG",
								"GGGTAAGTTGAGAGACAGGTTGGACAGGGTTAGATTAGGGCTGTGTTAGGGTAGTGTTAGGATGTGTGTGTGTGGGTGTGGTGTGGTGTGTGGTGTGGTG",
								"GGGTAAGTTGAGAGACAGGTTGGACAGGGTTAGATTAGGGCTGTGTTAGGGTAGTGTTAGGATGTGTGTGTGTGGGTGTGGTGTGGTGTGTGGTGTGGTG" };
		char const
	 * refs[] =
						{
								"GTTGGGTGACACACCCACACACCACACCACACACCAGACCCACACCCACAAACACACATCCTAAGACTGCCCTAAAACTGCCCTAATCTAACCCTGGCCAACCTGTCTCTGTGGTCA",
								"TCAACTTCCCACACACCACACCACACACCAGACCCACACCCACACACACACATCCTAAGACTGCCCTAAAACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCCCTCCAT",
								"AGGGTAACGCACACCCACACACCACACCACACACCAGACCCACACCCACACACACACATCCTAAGACTGCCCTAAAACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTTTGGGTG",
								"TTGGAGGGCACACCCACACACCACACCACACACCAGACCCACACCCACACACAAAACACATCCTAAGACTGCCCTAAAACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTTAACTTTG",
								"CCCTCCATTACACACCCACACACCACACCACACACCAGACCCACACCCAAGGACACACATCCTAAGACTGCCCTAAAACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCCCTGC",
								"TAATTGGAGGGTAAGTTGAGAGACAGGTTGGACAGGGTTAGATTAGGGCTGTGTTAGGGTAGTGTTAGGATGTGTGTGTGTGGGTGTGGTGTGGTGTGTGGTGTGGGTAACG",
								"TCCATTAAAGTTGAGAGACAGGTTGGACAGGGTTAGATTAGGGCTGTGTTAGGGTAGTGTTAGGATGTGTGTGTGTGGGTGTGGTGTGGTGTGTGGTGTGGTGCCCTGCCTG",
								"TATTACCCTGGGGTAAGTTGAGAGACAGGTTGGACAGGGTTAGATTAGGGCTGTGTTAGGGTAGTGTTAGGATGTGTGTGTGTGGGTGTGGTGTGGTGTGTGGTGTGGTGCCTCCA" };*/

//	char const
//	* refs[] =
//	{
//			"TTAATTTT",
//			"AAAAAAAA",
//			"AAAAAAAA",
//			"AAAAAAAA",
//			"AAAAAAAA",
//			"AAAAAAAA",
//			"AAAAAAAA",
//			"ATATTATA" };
//	char const
//	* reads[] =
//	{
//			"aa",
//			"AAAAAAAA",
//			"AATTTTAA",
//			"AAAATAAA",
//			"TTTTAAAA",
//			"AAAATTTT",
//			"ATATATAT",
//			"ATATATAT" };

//	char const
//	* refs[] =
//	{
//			"TTTAA",
//			"TTTGGCCTT",
//			"GGGGGTTTT",
//			"AACCCCCCAA",
//			"ATGC",
//			"AAAAAAAA",
//			"AAAAAAAA",
//			"ATATTATA" };
//	char const
//	* reads[] =
//	{
//			"aa",
//			"GGCC",
//			"TTTTAAA",
//			"AACCGGCCAA",
//			"ATGC",
//			"AAAATTTT",
//			"ATATATAT",
//			"ATATATAT" };

	char const
	* reads[] =
	{
			"AGGGGGGA",
			"AATTTTGCC",
			"TTTTTAA",
			"ATAGATAGATAG",
			"AGCAGTAC",
			"AGAGAGAG",
			"",
			"ATATATAT" };

	char const
	* refs[] =
	{
			"AGCAGATGAGGGCGGATAGC",
			"TTTTGCCAACGCATGGCAGA",
			"ATGACGACGCAGTGCTTTTT",
			"GCAATAGAATAGATAGTGGT",
			"TAGCATCAAGCAGTACTACA",
			"TCTCTCTCTCTCTCTCTCTC",
			"AGCAGATGACATGCATGCAA",
			"" };

	int seqNumber = 8;

	size_t max_read_length = pad(reads, seqNumber, READ_PAD);

	size_t max_ref_length =	pad(refs, seqNumber, REF_PAD);

	SSEKernel * ssekernel = new SSEKernel();

	short * scores = new short[seqNumber];

	ssekernel->init(max_read_length, max_ref_length);

	Timer timer;

	timer.start();

	Alignment * alignments = 0;

	for (int i = 0; i < 1; ++i) {

//		ssekernel->score_alignment_needleman_wunsch(reads, refs, scores);
//		for (int j = 0; j < seqNumber; ++j) {
//			cout << "Read:\t" << reads[j] << std::endl <<
//					"Ref:\t" << refs[j] << std::endl <<
//					"Score:\t" << scores[j] << std::endl;
//		}
		alignments = new Alignment[seqNumber];
		//ssekernel->calc_alignment(reads, refs, alignments);
		ssekernel->calc_alignment_needleman_wunsch(reads, refs, alignments);

	}

	timer.stop();

	for (int i = 0; i < 8; ++i) {
		std::cout << "==================" << std::endl << "\"";
		std::cout << alignments[i].read + alignments[i].readStart;
		std::cout << "\"" << std::endl << "\"";
		std::cout << alignments[i].ref + alignments[i].refStart;
		std::cout << "\"" << std::endl << "==================" << std::endl;
	}

	cout << "Alignment took " << timer.getElapsedTimeInMicroSec() / 10000 << " ms" << endl;

	// Premature return for SSE testing

//	for (int i = 0; i < seqNumber; ++i) {
//		cout << reads[i] << ":\t" << scores[i] << endl;
//	}

	delete ssekernel; ssekernel = 0;
	delete scores; scores = 0;

	SWKernel * kernel = new SWKernel();
	kernel->init(max_read_length, max_ref_length);

//	for (int i = 0; i < seqNumber; ++i) {
//
//		short * const score = (short * const)malloc(sizeof(short));
//		memset(score, 0, sizeof(short));
//
//		char const * const * const read = reads + i;
//		char const * const * const ref = refs + i;
//
//		timer.start();
//
//		for (int j = 0; j < 1; ++j) {
//
//			//for (int k = 0; k < 8; ++k) {
//
//				//kernel->score_alignment(read, ref, score);
//				kernel->score_alignment_needleman_wunsch(read, ref, score);
//			//}
//
//		}
//
//		cout << *(read) << std::endl << *(ref) << std::endl << "Score:\t" << *score << endl;
//		free(score);
//
//		timer.stop();
//
//		cout << "Alignment took " << timer.getElapsedTimeInMicroSec() / 10000 << " ms" << endl;
//	}

	alignments = new Alignment[seqNumber];

	for (int i = 0; i < seqNumber; ++i) {
		char const * const * const read = reads + i;
		char const * const * const ref = refs + i;

		std::cout << "Read: " << *read << std::endl;
		std::cout << "Ref: " << *ref << std::endl;

		timer.start();

		kernel->calc_alignment_needleman_wunsch(read, ref, &alignments[i]);

		timer.stop();

		cout << "Alignment took " << timer.getElapsedTimeInMicroSec() << " ms" << endl;

		std::cout << "==================" << std::endl << "\"";
		std::cout << alignments[i].read + alignments[i].readStart;
		std::cout << "\"" << std::endl << "\"";
		std::cout << alignments[i].ref + alignments[i].refStart;
		std::cout << "\"" << std::endl << "==================" << std::endl;

	}

	delete kernel; kernel = 0;

	cout << "Sizeof char * " << sizeof(char *) << endl;
	cout << "Sizeof int * " << sizeof(int *) << endl;
	cout << "Sizeof float " << sizeof(float) << endl;
	cout << "Sizeof int " << sizeof(int) << endl;
	cout << "Sizeof short " << sizeof(short) << endl;

	char N = 78;
	char n = 110;

	char mask = 223;

	char masked = n & mask;

	cout << N << " " << n << " " << mask << " " << masked;
	// N 01001110 -> 78
	// n 01101110 -> 110

	// Mask 11011111 -> 223

	return 0;
}

