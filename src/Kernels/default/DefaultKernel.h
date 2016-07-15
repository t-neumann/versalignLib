/*
 * DefaultKernel.h
 *
 *  Created on: May 23, 2016
 *      Author: tobias.neumann
 */

#ifndef DEFAULTKERNEL_H
#define DEFAULTKERNEL_H

#include "AlignmentKernel.h"
#include "AlignmentParameters.h"

#include <cstring>
#include <cmath>
#include <iostream>
#include <xmmintrin.h>
#include <stdlib.h>
#include <cstring>
#include <ctime>
#include <climits>

#define ASCII_ALPHABET 256
#define SCORE_CASE 6

#define UP 'u'
#define LEFT 'l'
#define DIAG 'd'
#define START 's'

typedef char * alnMat;

// Match anything non AGCTNagctn to 0
// Aa = 1
// Tt = 2
// Cc = 3
// Gg = 4
// Nn = 5

const char char_to_score[ASCII_ALPHABET] = {
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 1, 0, 3, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 5, 0,
		0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 1, 0, 3, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 5, 0,
		0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

class DefaultKernel: public AlignmentKernel {

public:

	DefaultKernel() {

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

//		std::cout << "Match: " << scoreMatch
//				<< "\nMismatch: " << scoreMismatch
//				<< "\nGap_read: " << scoreGapRead
//				<< "\nGap_ref: " << scoreGapRef
//				<< "\nRead_length: " << readLength
//				<< "\nRef_length: " << refLength
//				<< "\nAln_length: " << alnLength << std::endl;

		short tmp[SCORE_CASE][SCORE_CASE]= {
				// non ATGCN
				{0,0,0,0,0,0},
				// A
				{0,scoreMatch,scoreMismatch,scoreMismatch,scoreMismatch,0},
				// T
				{0,scoreMismatch,scoreMatch,scoreMismatch,scoreMismatch,0},
				// C
				{0,scoreMismatch,scoreMismatch,scoreMatch,scoreMismatch,0},
				// G
				{0,scoreMismatch,scoreMismatch,scoreMismatch,scoreMatch,0},
				// N
				{0,0,0,0,0,0}
		};
		memcpy(base_score, tmp, SCORE_CASE * SCORE_CASE * sizeof(short));


	}

	virtual ~DefaultKernel() {}

	virtual void compute_alignments(int const & opt, int const & aln_number, char const * const * const reads,
			char const * const * const refs, Alignment * const alignments);

	virtual void score_alignments(int const & opt, int const & aln_number, char const * const * const reads,
			char const * const * const refs, short * const scores);

private:

	virtual void score_alignment_smith_waterman(char const * const read,
			char const * const ref, short * const scores);

	virtual void score_alignment_needleman_wunsch(char const * const read,
			char const * const ref, short * const scores);

	virtual void calc_alignment_smith_waterman(char const * const read,
			char const * const ref, Alignment * const alignment);

	virtual void calc_alignment_needleman_wunsch(char const * const read,
			char const * const ref, Alignment * const alignment);

	void calculate_alignment_matrix_smith_waterman(char const * const read,
			char const * const ref, alnMat const matrix, short * const best_coordinates);

	void calculate_alignment_matrix_needleman_wunsch(char const * const read,
			char const * const ref, alnMat const matrix, short * const best_coordinates);

	int readLength;
	int refLength;
	int alnLength;

	short base_score[SCORE_CASE][SCORE_CASE];

	short scoreMatch;
	short scoreMismatch;
	short scoreGapRead;
	short scoreGapRef;

	typedef void (DefaultKernel::* fp_alignment_call)(char const * const,  char const * const, Alignment * const);
	typedef void (DefaultKernel::* fp_scoring_call)(char const * const,  char const * const, short * const);

};

#endif /* DEFAULTKERNEL_H */
