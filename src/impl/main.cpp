/*
 * main.cpp
 *
 *  Created on: May 20, 2016
 *      Author: tobias.neumann
 */

#include "AlignmentKernel.h"
#include "CustomParameters.h"
#include "dll_handling.h"
//#include "Kernels/plain/SWKernel.h"
#include "Kernels/AVX-SSE/SSEKernel.h"
//#include "Kernels/AVX-SSE/AVXKernel.h"
#include "util/versalignUtil.h"
#include "Kernels/OpenCL/ocl_testing.h"
#include "Timer/Timer.h"
#include <iostream>

#define READ_PAD '\0'
#define REF_PAD '\0'

using std::cout;
using std::string;
using std::endl;

CustomParameters parameters;

int main(int argc, char *argv[]) {

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
			"ATATATAT",
			"AGG",
			"TAA",
			"AGAG",
			"AGCAGTAG",
			"TATAC",
			"AGGAAGAG",
			"TT",
			"CCCCC"};

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
			"",
			"AGCAGATGAGGGCGGATAGC",
			"TTTTGCCAACGCATGGCAGA",
			"ATGACGACGCAGTGCTTTTT",
			"GCAATAGAATAGATAGTGGT",
			"TAGCATCAAGCAGTACTACA",
			"TCTCTCTCTCTCTCTCTCTC",
			"AGCAGATGACATGCATGCAA",
			""};

	int seqNumber = 16;

	size_t max_read_length = pad(reads, seqNumber, READ_PAD);

	size_t max_ref_length =	pad(refs, seqNumber, REF_PAD);

	parameters.read_length = max_read_length;
	parameters.ref_length = max_ref_length;

	char const * libPath = "../bin/libDefaultKernel.so";

	int const dll = DLL_init(libPath, &parameters);

	fp_load_alignment_kernel load_alignment_kernel = (fp_load_alignment_kernel) DLL_function_retreival(dll, "spawn_alignment_kernel");
	fp_delete_alignment_kernel delete_alignment_kernel = (fp_delete_alignment_kernel) DLL_function_retreival(dll, "delete_alignment_kernel");

	AlignmentKernel * plain_kernel = 0;

	plain_kernel = load_alignment_kernel();

	short * scores = new short[seqNumber]();
	Alignment * alignments = new Alignment[seqNumber]();

	plain_kernel->compute_alignments(0, seqNumber, reads, refs, alignments);
	plain_kernel->score_alignments(0, seqNumber, reads,refs, scores);

	//SWKernel * kernel = new SWKernel();

	//kernel->compute_alignments(0,seqNumber,reads,refs,alignments);

		for (int i = 0; i < seqNumber; ++i) {
			std::cout << "Read: " << reads[i] << std::endl;
			std::cout << "Ref: " << refs[i] << std::endl;
			std::cout << "==================" << std::endl << "\"";
			std::cout << alignments[i].read + alignments[i].readStart;
			std::cout << "\"" << std::endl << "\"";
			std::cout << alignments[i].ref + alignments[i].refStart;
			std::cout << "\"" << std::endl << "==================" << std::endl;
			std::cout << "Score: " << scores[i] << std::endl << std::endl;

		}

//	run_ocl_test(reads, refs, seqNumber, max_read_length, max_ref_length);

//	return 0;

//	AVXKernel * avxkernel = new AVXKernel();
//
//	short * avxscores = new short[seqNumber];
//
//	avxkernel->init(max_read_length, max_ref_length);
//
//	Alignment * avx_alignment = new Alignment[seqNumber];
////	avxkernel->score_alignment_needleman_wunsch(reads, refs, avxscores);
////
////	for (int j = 0; j < seqNumber; ++j) {
////		cout << "Read:\t" << reads[j] << std::endl <<
////				"Ref:\t" << refs[j] << std::endl <<
////				"Score:\t" << avxscores[j] << std::endl;
////	}
////	std::cout << std::endl << "AVX end" << std::endl;
//
//	avxkernel->calc_alignment_needleman_wunsch(reads, refs, avx_alignment);
//
//	for (int j = 0; j < seqNumber; ++j) {
//		std::cout << "==================" << std::endl << "\"";
//		std::cout << avx_alignment[j].read + avx_alignment[j].readStart;
//		std::cout << "\"" << std::endl << "\"";
//		std::cout << avx_alignment[j].ref + avx_alignment[j].refStart;
//		std::cout << "\"" << std::endl << "==================" << std::endl;
//	}
//
//	delete avxkernel;

//	SSEKernel * ssekernel = new SSEKernel();
//
//	short * scores = new short[seqNumber];
//
//	ssekernel->init(max_read_length, max_ref_length);
//
//	Timer timer;
//
//	timer.start();
//
//	Alignment * alignments = 0;
//
//	for (int i = 0; i < 2; ++i) {
//
//		char * * reads_batch = new char * [8];
//		for (int j = 0; j < 8; ++j) {
//			reads_batch[j] = new char [max_read_length];
//			memcpy(reads_batch[j],reads[i * 8 + j], max_read_length * sizeof(char));
//			//reads_batch[j] = reads[i * 8 + j];
//		}
//		char * * refs_batch = new char * [8];
//		for (int j = 0; j < 8; ++j) {
//			refs_batch[j] = new char [max_ref_length];
//			memcpy(refs_batch[j],refs[i * 8 + j], max_ref_length * sizeof(char));
//			//refs_batch[j] = refs[i * 8 + j];
//		}
//
//		//ssekernel->score_alignment_needleman_wunsch(reads, refs, scores);
//		ssekernel->score_alignment_needleman_wunsch(reads_batch, refs_batch, scores);
//		for (int j = 0; j < 8; ++j) {
//			cout << "Read:\t" << reads_batch[j] << std::endl <<
//					"Ref:\t" << refs_batch[j] << std::endl <<
//					"Score:\t" << scores[j] << std::endl;
//		}
////		std::cout << std::endl << "Batch end" << std::endl;
////		alignments = new Alignment[8];
////		ssekernel->calc_alignment_needleman_wunsch(reads_batch, refs_batch, alignments);
//		//ssekernel->calc_alignment_needleman_wunsch(reads, refs, alignments);
////		ssekernel->calc_alignment_needleman_wunsch(reads_batch, refs_batch, alignments);
//
////		for (int j = 0; j < 8; ++j) {
////			std::cout << "==================" << std::endl << "\"";
////			std::cout << alignments[j].read + alignments[j].readStart;
////			std::cout << "\"" << std::endl << "\"";
////			std::cout << alignments[j].ref + alignments[j].refStart;
////			std::cout << "\"" << std::endl << "==================" << std::endl;
////		}
////
//	}
//
//	timer.stop();
//
//	cout << "Alignment took " << timer.getElapsedTimeInMicroSec() / 10000 << " ms" << endl;
//
//	return 0;

	// Premature return for SSE testing

//	for (int i = 0; i < seqNumber; ++i) {
//		cout << reads[i] << ":\t" << scores[i] << endl;
//	}

//	delete ssekernel; ssekernel = 0;
//	delete scores; scores = 0;

//	SWKernel * kernel = new SWKernel();
	//kernel->init(max_read_length, max_ref_length);

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

//	alignments = new Alignment[seqNumber];
//
//	for (int i = 0; i < seqNumber; ++i) {
//		char const * const * const read = reads + i;
//		char const * const * const ref = refs + i;
//
//		std::cout << "Read: " << *read << std::endl;
//		std::cout << "Ref: " << *ref << std::endl;
//
//		timer.start();
//
//		kernel->calc_alignment_needleman_wunsch(read, ref, &alignments[i]);
//
//		timer.stop();
//
//		cout << "Alignment took " << timer.getElapsedTimeInMicroSec() << " ms" << endl;
//
//		std::cout << "==================" << std::endl << "\"";
//		std::cout << alignments[i].read + alignments[i].readStart;
//		std::cout << "\"" << std::endl << "\"";
//		std::cout << alignments[i].ref + alignments[i].refStart;
//		std::cout << "\"" << std::endl << "==================" << std::endl;
//
//	}

//	delete kernel; kernel = 0;

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

