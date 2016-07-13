/*
 * AVXKernel.h
 *
 *  Created on: May 23, 2016
 *      Author: tobias.neumann
 */

#ifndef AVXKERNEL_H_
#define AVXKERNEL_H_

//#include "AlignmentKernel.h"
#include "AlignmentKernel.h"

#if _WIN32
#define align32 __declspec(align(32))
#define malloc32(ptr,size,align)  ptr = (short*) _aligned_malloc(size, align)
#else
#include <stdlib.h>
#define malloc32(ptr,size,align)  posix_memalign(((void * *)&ptr), align, size)
#define align32 __attribute__((aligned(32)))
#endif

#define AVX_SIZE 16

#define UP 1
#define LEFT 2
#define DIAG 3
#define START 0

#include <immintrin.h>
#include <climits>

class AVXKernel {

public:
	AVXKernel() {

		scoreGapRead = -3;
		scoreGapRef = -3;
		scoreMatch = 2;
		scoreMismatch = -1;

		x_scoreMatch = short_to_avx(scoreMatch);
		x_scoreMismatch = short_to_avx(scoreMismatch);
		x_scoreGapRead = short_to_avx(scoreGapRead);
		x_scoreGapRef = short_to_avx(scoreGapRef);

		x_zeros = short_to_avx(0);

		x_A = short_to_avx(65);
		x_T = short_to_avx(84);
		x_C = short_to_avx(67);
		x_G = short_to_avx(71);

		x_UCMask = short_to_avx(223);

		x_p_up = short_to_avx(UP);
		x_p_left = short_to_avx(LEFT);
		x_p_diag = short_to_avx(DIAG);
	}

	void init (int const & max_read_length, int const & max_ref_length) {
			this->readLength = max_read_length;
			this->refLength = max_ref_length;
			this->alnLength = refLength + readLength;
	}

	virtual ~AVXKernel() {}

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
	// __m256i fits 256 bits = 16 shorts
	inline __m256i short_to_avx(short x) {

		align32 short buf[AVX_SIZE];
		for (int i = 0; i < AVX_SIZE; ++i)
			buf[i] = x;

		return _mm256_load_si256((__m256i *) buf);
	}

	inline __m256i _mm_blendv_si256 (__m256i x, __m256i y, __m256i mask) {
		// Replace bit in x with bit in y when matching bit in mask is set:
		return _mm256_or_si256(_mm256_andnot_si256(mask, x), _mm256_and_si256(mask, y));
	}

	int readLength;
	int refLength;
	int alnLength;

	short scoreMatch;
	short scoreMismatch;
	short scoreGapRead;
	short scoreGapRef;

	__m256i x_scoreMatch;
	__m256i x_scoreMismatch;
	__m256i x_scoreGapRead;
	__m256i x_scoreGapRef;

	__m256i x_zeros;

	__m256i x_A;


	// A 01000001 -> 65
	// a 01100001 -> 97

	__m256i x_T;

	// T 01010100 -> 84
	// t 01110100 -> 116

	__m256i x_C;

	// C 01000011 -> 67
	// c 01100011 -> 99

	__m256i x_G;

	// G 01000111 -> 71
	// g 01100111 -> 103

	// N 01001110 -> 78
	// n 01101110 -> 110

	// Mask 11011111 -> 223

	__m256i x_UCMask;

	// Pointer directions
	__m256i x_p_up;
	__m256i x_p_left;
	__m256i x_p_diag;

};

#endif /* AVXKERNEL_H_ */
