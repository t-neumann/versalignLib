/*
 * SSEKernel.h
 *
 *  Created on: May 23, 2016
 *      Author: tobias.neumann
 */

#ifndef SSEKERNEL_H_
#define SSEKERNEL_H_

class SSEKernel: public AlignmentKernel {
public:
	SSEKernel();
	virtual ~SSEKernel();

	float score_alignment(char const * read, float const & read_length,
			char const * reference, float const & reference_length,
			float const & gap_read, float const & gap_ref, float const & match,
			float const & mismatch) = 0;
	float score_alignment_corridor(char const * read, float const & read_length,
			char const * reference, float const & reference_length,
			float const & gap_read, float const & gap_ref, float const & match,
			float const & mismatch, float const & corridor_width) = 0;
};

#endif /* SSEKERNEL_H_ */
