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

	free(tmp);
}

void SSEKernel::score_alignment (char const * const * const read, char const * const * const ref, short * const scores) {

	__m128i max_score = short_to_sse(0);

	// Initialize SSE matrix
	short * matrix = 0;

	malloc16(matrix, sizeof(short) * (refLength + 1) * SCORING_ROWS * SSE_SIZE,16);
	memset(matrix, 0, sizeof(short) * (refLength + 1) * SCORING_ROWS * SSE_SIZE);

	// offset for first and second row
	int prev_row = 0;
	int cur_row = 1;

	// SSE conversion array for current read and ref bases
	align16 short read_bases [SSE_SIZE];
	align16 short ref_bases  [SSE_SIZE];

	// holds current read and ref base for SSE instruction
	__m128i sse_read_bases;
	__m128i sse_ref_bases;

	for (int read_pos = 0; read_pos < readLength; ++read_pos) {

		// load read base
		for (int base = 0; base < SSE_SIZE; ++base) {
			read_bases[base] = *(read[base] + read_pos);
		}

		sse_read_bases = _mm_load_si128((__m128i const *) read_bases);

//		std::cout << "Read base:\t";
//		print_sse_char(sse_read_bases);

		for (int ref_pos = 0; ref_pos < refLength; ++ref_pos) {

			// load ref base
			for (int base = 0; base < SSE_SIZE; ++base) {
				ref_bases[base] = *(ref[base] + ref_pos);
			}

			sse_ref_bases = _mm_load_si128((__m128i const *) ref_bases);

//			std::cout << "Ref base:\t";
//			print_sse_char(sse_ref_bases);

			// load relevant matrix cells (up, diag, left)
			__m128i up = _mm_load_si128((__m128i *) (matrix + SSE_SIZE * (prev_row * (refLength + 1) + ref_pos + 1)));
			__m128i diag = _mm_load_si128((__m128i *) (matrix + SSE_SIZE * (prev_row * (refLength + 1) + ref_pos)));
			__m128i left = _mm_load_si128((__m128i *) (matrix + SSE_SIZE * (cur_row * (refLength + 1) + ref_pos)));

//			std::cout << "Up:\t";
//			print_sse_short(up);
//			std::cout << "Left:\t";
//			print_sse_short(left);
//			std::cout << "Diag:\t";
//			print_sse_short(diag);

			// add gap penalties to up and left
			up = _mm_add_epi16(up, x_scoreGapRef);
			left = _mm_add_epi16(left, x_scoreGapRead);

			// match read and ref bases
			// matches will hold 1, mismatches will hold 0
			__m128i match = _mm_cmpeq_epi16(sse_read_bases, sse_ref_bases);

			// bitwise and between 0 (mismatches) and match score gives match score only for equal bases
			// add these to diagonal value
			diag = _mm_add_epi16(diag, _mm_and_si128(match, x_scoreMatch));
			// bitwise not and between 1 (matches) and mismatch score gives mismatch score only for unequal bases
			// add these to diagonal value
			diag = _mm_add_epi16(diag, _mm_andnot_si128(match, x_scoreMismatch));

			// Cell value will be max of upper + gap penalty, left + gap penalty, diag + match/mismatch score or 0
			__m128i cell = _mm_max_epi16(diag, _mm_max_epi16(left, _mm_max_epi16(up, x_zeros)));

//			std::cout << "Up:\t";
//			print_sse_short(up);
//			std::cout << "Left:\t";
//			print_sse_short(left);
//			std::cout << "Diag:\t";
//			print_sse_short(diag);
//			std::cout << "Max:\t";
//			print_sse_short(cell);

			_mm_store_si128((__m128i *) (matrix + SSE_SIZE * (cur_row * (refLength + 1) + ref_pos + 1)), cell);

			max_score = _mm_max_epi16(cell, max_score);

		}
		prev_row = cur_row;
		(++cur_row) &= 1;
	}

	free(matrix);

	_mm_store_si128((__m128i *)scores, max_score);
}

#undef SCORING_ROWS
#undef SSE_SIZE
