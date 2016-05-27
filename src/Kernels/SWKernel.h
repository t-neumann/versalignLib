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

struct mat_element;
struct xy_coordinates;

enum aln_paths {
	start,
	upwards,
	leftwards,
	diagonal
};

class SWKernel: public AlignmentKernel {

public:

	SWKernel() {

		scoreGapRead = -3;
		scoreGapRef = -3;
		scoreMatch = 2;
		scoreMismatch = -1;

	}

	virtual ~SWKernel();

	void set_reference_length(int const & reference_size) {
		this->refLength = reference_size;
	}

	void set_read_length(int const & read_length) {
		this->readLength = read_length;
	}

	virtual void score_alignment(char const * const * const read,
			char const * const * const ref, short * const scores);

	virtual void calc_alignment(char const * const * const read,
			char const * const * const ref, char * const * const aligned_read,
			char * const * const aligned_ref);

private:

	xy_coordinates calculate_alignment_matrix(char const * const * const read,
			char const * const * const ref, mat_element * const matrix);
	int readLength;
	int refLength;

	short scoreMatch;
	short scoreMismatch;
	short scoreGapRead;
	short scoreGapRef;
};

struct xy_coordinates {
	short x;
	short y;
};

struct mat_element {
	short score;
	aln_paths path;

	mat_element() {
		score = 0;
		path = start;
	}
};

#endif /* SWKERNEL_H_ */
