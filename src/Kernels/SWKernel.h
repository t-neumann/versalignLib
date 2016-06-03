/*
 * SWKernel.h
 *
 *  Created on: May 23, 2016
 *      Author: tobias.neumann
 */

#ifndef SWKERNEL_H_
#define SWKERNEL_H_

#include "AlignmentKernel.h"

#include <cstring>
#include <cmath>
#include <iostream>
#include <xmmintrin.h>
#include <stdlib.h>
#include <cstring>
#include <ctime>

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

class SWKernel: public AlignmentKernel {

public:

	SWKernel() {

		scoreGapRead = -3;
		scoreGapRef = -3;
		scoreMatch = 2;
		scoreMismatch = -1;

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

	virtual ~SWKernel();

	void init (int const & max_read_length, int const & max_ref_length) {
		this->readLength = max_read_length;
		this->refLength = max_ref_length;
		this->alnLength = refLength + readLength;
	}

	void set_reference_length(int const & reference_size) {
		this->refLength = reference_size;
	}

	void set_read_length(int const & read_length) {
		this->readLength = read_length;
	}

	virtual void score_alignment(char const * const * const read,
			char const * const * const ref, short * const scores);

	virtual void calc_alignment(char const * const * const read,
			char const * const * const ref, Alignment * const alignment);

private:

	void calculate_alignment_matrix(char const * const * const read,
			char const * const * const ref, alnMat const matrix, short * const best_coordinates);
	int readLength;
	int refLength;
	int alnLength;

	short base_score[SCORE_CASE][SCORE_CASE];

	short scoreMatch;
	short scoreMismatch;
	short scoreGapRead;
	short scoreGapRef;
};

#endif /* SWKERNEL_H_ */
