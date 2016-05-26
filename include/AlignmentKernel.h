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

	virtual float score_alignment(char const * read, float const & read_length,
			char const * reference, float const & reference_length,
			float const & gap_read, float const & gap_ref, float const & match,
			float const & mismatch) = 0;
	virtual float score_alignment_corridor(char const * read,
			float const & read_length, char const * reference,
			float const & reference_length, float const & gap_read,
			float const & gap_ref, float const & match, float const & mismatch,
			float const & corridor_width) = 0;
};

#endif /* INCLUDE_ALIGNMENTKERNEL_H_ */
