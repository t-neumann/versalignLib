/*
 * SSEKernel.h
 *
 *  Created on: May 23, 2016
 *      Author: tobias.neumann
 */

#ifndef SSEKERNEL_H_
#define SSEKERNEL_H_

//#include "AlignmentKernel.h"
#include "../../include/AlignmentKernel.h"

#if _WIN32
#define align16 __declspec(align(16))
#define malloc16(ptr,size,align)  ptr = (short*) _aligned_malloc(size, align)
#else
#include <stdlib.h>
#define malloc16(ptr,size,align)  posix_memalign(((void * *)&ptr), align, size)
#define align16 __attribute__((aligned(16)))
#endif

#include <emmintrin.h>

class SSEKernel: public AlignmentKernel {

public:
	SSEKernel() {

		scoreGapRead = -3;
		scoreGapRef = -3;
		scoreMatch = 2;
		scoreMismatch = -1;

		x_scoreMatch = short_to_sse(scoreMatch);
		x_scoreMismatch = short_to_sse(scoreMismatch);
		x_scoreGapRead = short_to_sse(scoreGapRead);
		x_scoreGapRef = short_to_sse(scoreGapRef);

		x_zeros = short_to_sse(0);
	}

	virtual ~SSEKernel() {}

	void set_reference_length(int const & reference_size) {
		this->refLength = reference_size;
	}

	void set_read_length(int const & read_length) {
		this->readLength = read_length;
	}

	void score_alignment(char const * const * const read,
			char const * const * const ref, short * const scores);

private:

	// Short = 2 byte
	// __m128i fits 128 bits = 8 shorts

	inline __m128i short_to_sse(short x) {

		align16 short buf[8];
		for (int i = 0; i < 8; ++i)
			buf[i] = x;

		return _mm_load_si128((__m128i *) buf);
	}

	int readLength;
	int refLength;

	short scoreMatch;
	short scoreMismatch;
	short scoreGapRead;
	short scoreGapRef;

	__m128i x_scoreMatch;
	__m128i x_scoreMismatch;
	__m128i x_scoreGapRead;
	__m128i x_scoreGapRef;

	__m128i x_zeros;

};

#endif /* SSEKERNEL_H_ */
