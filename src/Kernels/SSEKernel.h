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

#define SSE_SIZE 8

#define UP 1
#define LEFT 2
#define DIAG 3
#define START 0

#include <emmintrin.h>
#include <climits>

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

		x_A = short_to_sse(65);
		x_T = short_to_sse(84);
		x_C = short_to_sse(67);
		x_G = short_to_sse(71);

		x_UCMask = short_to_sse(223);

		x_p_up = short_to_sse(UP);
		x_p_left = short_to_sse(LEFT);
		x_p_diag = short_to_sse(DIAG);
	}

	void init (int const & max_read_length, int const & max_ref_length) {
			this->readLength = max_read_length;
			this->refLength = max_ref_length;
			this->alnLength = refLength + readLength;
	}

	virtual ~SSEKernel() {}

	void score_alignment(char const * const * const read,
			char const * const * const ref, short * const scores);

	void score_alignment_needleman_wunsch(char const * const * const read,
				char const * const * const ref, short * const scores);

	void calc_alignment(char const * const * const read,
				char const * const * const ref, Alignment * const alignment);

	virtual void calc_alignment_needleman_wunsch(char const * const * const read,
				char const * const * const ref, Alignment * const alignment);

private:

	void calc_alignment_matrix(char const * const * const read,
			char const * const * const ref, short * const matrix, short * const best_coordinates);

	void calculate_alignment_matrix_needleman_wunsch(char const * const * const read,
			char const * const * const ref, short * const matrix, short * const best_coordinates);

	// Short = 2 byte
	// __m128i fits 128 bits = 8 shorts
	inline __m128i short_to_sse(short x) {

		align16 short buf[SSE_SIZE];
		for (int i = 0; i < SSE_SIZE; ++i)
			buf[i] = x;

		return _mm_load_si128((__m128i *) buf);
	}

	inline __m128i _mm_blendv_si128 (__m128i x, __m128i y, __m128i mask) {
		// Replace bit in x with bit in y when matching bit in mask is set:
		return _mm_or_si128(_mm_andnot_si128(mask, x), _mm_and_si128(mask, y));
	}

	int readLength;
	int refLength;
	int alnLength;

	short scoreMatch;
	short scoreMismatch;
	short scoreGapRead;
	short scoreGapRef;

	__m128i x_scoreMatch;
	__m128i x_scoreMismatch;
	__m128i x_scoreGapRead;
	__m128i x_scoreGapRef;

	__m128i x_zeros;

	__m128i x_A;


	// A 01000001 -> 65
	// a 01100001 -> 97

	__m128i x_T;

	// T 01010100 -> 84
	// t 01110100 -> 116

	__m128i x_C;

	// C 01000011 -> 67
	// c 01100011 -> 99

	__m128i x_G;

	// G 01000111 -> 71
	// g 01100111 -> 103

	// N 01001110 -> 78
	// n 01101110 -> 110

	// Mask 11011111 -> 223

	__m128i x_UCMask;

	// Pointer directions
	__m128i x_p_up;
	__m128i x_p_left;
	__m128i x_p_diag;

};

#endif /* SSEKERNEL_H_ */
