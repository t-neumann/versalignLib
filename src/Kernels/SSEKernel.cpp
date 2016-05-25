/*
 * SSEKernel.cpp
 *
 *  Created on: May 23, 2016
 *      Author: tobias.neumann
 */

#include "SSEKernel.h"
#include <memory.h>
#include<iostream>
#include<cstdio>

#define SCORING_ROWS 2
#define SSE_SIZE 8

float SSEKernel::score_alignment(char const * read, float const & read_length,
		char const * reference, float const & reference_length,
		float const & gap_read, float const & gap_ref, float const & match,
		float const & mismatch) {
	return 0.f;
}
float SSEKernel::score_alignment_corridor(char const * read,
		float const & read_length, char const * reference,
		float const & reference_length, float const & gap_read,
		float const & gap_ref, float const & match, float const & mismatch,
		float const & corridor_width) {

	return 0.f;

}

void print_sse_char (__m128i sse) {
	short * tmp;
	malloc16(tmp, sizeof(short) * 8, 16);

	_mm_store_si128((__m128i *)tmp, sse);

	for (int i = 0; i < SSE_SIZE; ++i) {
		int character = tmp[i];
		std::cout << i << ":" << (char)character << " ";
	}
	std::cout << std::endl;
}

void print_sse_short (__m128i sse) {
	short * tmp;
	malloc16(tmp, sizeof(short) * 8, 16);

	_mm_store_si128((__m128i *)tmp, sse);

	for (int i = 0; i < SSE_SIZE; ++i) {
		int character = tmp[i];
		std::cout << i << ":" << tmp[i] << " ";
	}
	std::cout << std::endl;
}

void SSEKernel::score_alignment (char const * const * const ref, char const * const * const read, short * scores) {

	malloc16(scores, sizeof(short) * (refLength + 1) * SCORING_ROWS * SSE_SIZE,16);
	memset(scores, 0, sizeof(short) * (refLength + 1) * SCORING_ROWS * SSE_SIZE);

	int prev_row = 0;
	int cur_row = 1;

	align16 short read_bases [SSE_SIZE];
	align16 short ref_bases  [SSE_SIZE];

	__m128i sse_read_bases;
	__m128i sse_ref_bases;

	for (int read_pos = 0; read_pos < readLength; ++ read_pos) {

		for (int base = 0; base < SSE_SIZE; ++base) {
			read_bases[base] = *(read[base] + read_pos);
		}

		sse_read_bases = _mm_load_si128((__m128i const *) read_bases);

		//print_sse_char(sse_read_bases);
		//print_sse_short(sse_read_bases);

		for (int ref_pos = 0; ref_pos < refLength; ++ref_pos) {

			for (int base = 0; base < SSE_SIZE; ++base) {
				ref_bases[base] = *(ref[base] + ref_pos);
			}

			sse_ref_bases = _mm_load_si128((__m128i const *) ref_bases);

			__m128i up = _mm_load_si128((__m128i *) (scores + SSE_SIZE * (prev_row * (refLength + 1) + ref_pos + 1)));
			__m128i diag = _mm_load_si128((__m128i *) (scores + SSE_SIZE * (prev_row * (refLength + 1) + ref_pos)));
			__m128i left = _mm_load_si128((__m128i *) (scores + SSE_SIZE * (cur_row * (refLength + 1) + ref_pos)));

			up = _mm_add_epi16(up, x_scoreGapRef);
			left = _mm_add_epi16(up, x_scoreGapRead);

//			__m128i cmp = _mm_cmpeq_epi16(q, r);
//			nw = _mm_add_epi16(nw, _mm_and_si128(cmp, m_xMatch));
//			nw = _mm_add_epi16(nw, _mm_andnot_si128(cmp, m_xMismatch));
//
//			__m128i cur = _mm_max_epi16(n, _mm_max_epi16(w, _mm_max_epi16(nw, m_xZero)));
//
//			score_max = _mm_max_epi16(cur, score_max);
//
//			_mm_store_si128((__m128i *) (scores + 8 * (cl * (m_Corridor + 2) + i)), cur);

		}
		prev_row = cur_row;
		(++cur_row) &= 1;
	}
}

#undef SCORING_ROWS
#undef SSE_SIZE
