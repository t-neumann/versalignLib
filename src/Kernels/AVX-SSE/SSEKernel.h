/*
 * SSEKernel.h
 *
 *  Created on: May 23, 2016
 *      Author: tobias.neumann
 */

#ifndef SSEKERNEL_H
#define SSEKERNEL_H

#include "AlignmentKernel.h"
#include "AlignmentParameters.h"
#include "AlignmentLogger.h"

#if _WIN32
#define align16 __declspec(align(16))
#define malloc16(ptr,size,align)  ptr = (short*) _aligned_malloc(size, align)
#else
#include <stdlib.h>
#define malloc16(ptr,size,align)  posix_memalign(((void * *)&ptr), align, size)
#define align16 __attribute__((aligned(16)))
#endif

#define SSE_SIZE 8

#define KERNEL "SSE"

#define UP 1
#define LEFT 2
#define DIAG 3
#define START 0

#include <immintrin.h>
#include <climits>

#include <iostream>

class SSEKernel : public AlignmentKernel {

public:
	SSEKernel() {

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

#ifndef NDEBUG
		Logger.log(0, KERNEL, "Successfully instantiated SSE Kernel.");
#endif

	}

	virtual ~SSEKernel() {}

	virtual void compute_alignments(int const & opt, int const & aln_number, char const * const * const reads,
			char const * const * const refs, Alignment * const alignments);

	virtual void score_alignments(int const & opt, int const & aln_number, char const * const * const reads,
			char const * const * const refs, short * const scores);

private:

	void score_alignment_smith_waterman(char const * const * const read,
			char const * const * const ref, short * const scores);

	void score_alignment_needleman_wunsch(char const * const * const read,
			char const * const * const ref, short * const scores);

	void calc_alignment_smith_waterman(char const * const * const read,
			char const * const * const ref, Alignment * const alignment);

	virtual void calc_alignment_needleman_wunsch(char const * const * const read,
			char const * const * const ref, Alignment * const alignment);

	void calc_alignment_matrix_smith_waterman(char const * const * const read,
			char const * const * const ref, short * const matrix, short * const best_coordinates);

	void calculate_alignment_matrix_needleman_wunsch(char const * const * const read,
			char const * const * const ref, short * const matrix, short * const best_coordinates);

	typedef void (SSEKernel::* fp_alignment_call)(char const * const * const,  char const * const * const, Alignment * const);
	typedef void (SSEKernel::* fp_scoring_call)(char const * const * const,  char const * const * const, short * const);

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

#endif /* SSEKERNEL_H */
