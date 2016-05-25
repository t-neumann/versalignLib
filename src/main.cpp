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

using std::cout;
using std::string;
using std::endl;

int main(int argc, char *argv[]) {

	int const RUNS = 1000;

	std::cout << "Startup\n";

	string b = "aaaaaaaaaabbbbbbbaaaaaaaaaa";
	const char * b_s = b.c_str();

	string a = "aaaaaaaaaaaaaaaaaaaa";
	const char *a_s = a.c_str();

	AlignmentKernel * kernel = new SWKernel();

	cout << "Scoring read:\t" << a << endl;
	cout << "Scoring ref:\t" << b << endl;

	Timer timer;

	timer.start();

	float score = kernel->score_alignment(a_s, (float)a.size(), b_s, (float)b.size(),3.0f, 3.0f, 2.0f, 5.0f);

	timer.stop();

	cout << "Alignment scored " << score << endl;
	cout << "Alignment took " << timer.getElapsedTimeInMicroSec() << " ms" << endl;

	delete kernel; kernel = 0;

	char const
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
								"TATTACCCTGGGGTAAGTTGAGAGACAGGTTGGACAGGGTTAGATTAGGGCTGTGTTAGGGTAGTGTTAGGATGTGTGTGTGTGGGTGTGGTGTGGTGTGTGGTGTGGTGCCTCCA" };
	int seqNumber = 8;

	size_t max_read_length = pad(reads, seqNumber);

	size_t max_ref_length =	pad(refs, seqNumber);

	SSEKernel * ssekernel = new SSEKernel();

	short * scores = 0;

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

	delete ssekernel; ssekernel = 0;

	cout << "Sizeof char * " << sizeof(char *) << endl;
	cout << "Sizeof int * " << sizeof(int *) << endl;
	cout << "Sizeof float " << sizeof(float) << endl;
	cout << "Sizeof int " << sizeof(int) << endl;
	cout << "Sizeof short " << sizeof(short) << endl;

	return 0;
}

