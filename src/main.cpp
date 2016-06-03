/*
 * main.cpp
 *
 *  Created on: May 20, 2016
 *      Author: tobias.neumann
 */

#include "Kernels/SWKernel.h"
#include "Kernels/SSEKernel.h"
#include "util/versalignUtil.h"
#include "Timer/Timer.h"
#include <iostream>
#include <string>

#define READ_PAD '\0'
#define REF_PAD '$'

using std::cout;
using std::string;
using std::endl;

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

	char const
					* refs[] =
							{
									"AA",
									"AAAAAAAA",
									"AAAAAAAA",
									"AAAAAAAA",
									"AAAAAAAA",
									"AAAAAAAA",
									"AAAAAAAA",
									"ATATTATA" };
			char const
					* reads[] =
							{
									"aa",
									"AAAAAAAA",
									"AATTTTAA",
									"AAAATAAA",
									"TTTTAAAA",
									"AAAATTTT",
									"ATATATAT",
									"ATATATAT" };

	int seqNumber = 8;

	size_t max_read_length = pad(reads, seqNumber, READ_PAD);

	size_t max_ref_length =	pad(refs, seqNumber, REF_PAD);

	SSEKernel * ssekernel = new SSEKernel();

	short * scores = new short[seqNumber];

	ssekernel->set_read_length(max_read_length);
	ssekernel->set_reference_length(max_ref_length);
	ssekernel->score_alignment(refs, reads, scores);

//		for (int j = 0; j < maxReadLen; ++j) {
//			for (int i = 0; i < seqNumber; ++i) {
//				cout << *(*(reads + i) + j) << std::endl;
//			}
//		}
//
//		cout << "MaxLen: " << maxReadLen << ", " << maxRefLen << endl;

	for (int i = 0; i < seqNumber; ++i) {
		cout << reads[i] << ":\t" << scores[i] << endl;
	}

	delete ssekernel; ssekernel = 0;
	delete scores; scores = 0;

	SWKernel * kernel = new SWKernel();
	kernel->init(max_read_length, max_ref_length);

	Timer timer;

	timer.start();

//	Alignment * alignments = new Alignment[seqNumber];
//
//	for (int i = 0; i < seqNumber; ++i) {
//		char const * const * const read = reads + i;
//		char const * const * const ref = refs + i;
//
//		short * const score = (short * const)malloc(sizeof(short));
//		memset(score, 0, sizeof(short));
//
//		kernel->score_alignment(read, ref, score);
//
//		cout << *(read) << ":\t"
//				<< *score << endl;
//		free(score);
//
//		kernel->calc_alignment(read, ref, &alignments[i]);
//
//		std::cout << "==================" << std::endl << "\"";
//		std::cout << alignments[i].read + alignments[i].readStart;
//		std::cout << "\"" << std::endl << "\"";
//		std::cout << alignments[i].ref + alignments[i].refStart;
//		std::cout << "\"" << std::endl << "==================" << std::endl;
//
//	}

	timer.stop();

	cout << "Alignment took " << timer.getElapsedTimeInMicroSec() << " ms" << endl;

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

