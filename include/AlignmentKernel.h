/*
 * Aligner.h
 *
 *  Created on: May 23, 2016
 *      Author: tobias.neumann
 */

#ifndef INCLUDE_ALIGNMENTKERNEL_H_
#define INCLUDE_ALIGNMENTKERNEL_H_

class AlignmentKernel {

public:

	virtual ~AlignmentKernel() {}

	virtual void score_alignment(char const * const * const read,
			char const * const * const ref, short * const scores) = 0;
};

struct Alignment {
	char * read;
	char * ref;
	short readStart;
	short readEnd;
	short refStart;
	short refEnd;
};

#endif /* INCLUDE_ALIGNMENTKERNEL_H_ */
