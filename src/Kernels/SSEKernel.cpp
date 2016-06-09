/*
 * SSEKernel.cpp
 *
 *  Created on: May 23, 2016
 *      Author: tobias.neumann
 */

#include "SSEKernel.h"
#include <memory.h>
#include <iostream>
#include <cstdio>
#include <sstream>
#include "../util/versalignUtil.h"

#define SCORING_ROWS 2

void SSEKernel::calc_alignment_matrix(char const * const * const read,
		char const * const * const ref, short * const matrix, short * const best_coordinates) {

	// Tracking best read and ref positions
	__m128i best_read_pos = x_zeros;
	__m128i best_ref_pos = x_zeros;
	__m128i max_score = x_zeros;

	// Scoring matrix indices
	short prev_row_score = 0;
	short current_row_score = 1;

	// Alignment matrix index
	short current_row_aln = 1;

	short * scoreMat = 0;

	malloc16(scoreMat, sizeof(short) * (refLength + 1) * SCORING_ROWS * SSE_SIZE,16);
	memset(scoreMat, 0, sizeof(short) * (refLength + 1) * SCORING_ROWS * SSE_SIZE);

	// SSE conversion array for current read and ref bases
	align16 short read_bases [SSE_SIZE];
	align16 short ref_bases  [SSE_SIZE];

	// holds current read and ref base for SSE instruction
	__m128i sse_read_bases;
	__m128i sse_ref_bases;

	// holds current read and ref positions
	__m128i sse_read_pos = x_zeros;
	__m128i sse_ref_pos = x_zeros;

	// holds 1 for increment
	__m128i see_increment = short_to_sse(1);

	for (int read_pos = 0; read_pos < readLength; ++read_pos) {

		sse_ref_pos = x_zeros;

		// load read base
		for (int base = 0; base < SSE_SIZE; ++base) {
			read_bases[base] = *(read[base] + read_pos);
		}

		sse_read_bases = _mm_load_si128((__m128i const *) read_bases);

		// UC read
		sse_read_bases = _mm_and_si128(sse_read_bases,x_UCMask);

		//std::cout << "Read base current:\t" << __m128i_toString<char>(sse_read_bases);

		__m128i valid_read_base = _mm_or_si128(_mm_cmpeq_epi16(sse_read_bases,x_C),_mm_or_si128(_mm_cmpeq_epi16(sse_read_bases,x_G),_mm_or_si128(_mm_cmpeq_epi16(sse_read_bases,x_T),_mm_cmpeq_epi16(sse_read_bases,x_A))));

		//std::cout << std::endl << "Valid:\t" << __m128i_toString<short>(valid_read_base) << std::endl;

		for (int ref_pos = 0; ref_pos < refLength; ++ref_pos) {

			// load ref base
			for (int base = 0; base < SSE_SIZE; ++base) {
				ref_bases[base] = *(ref[base] + ref_pos);
			}

			sse_ref_bases = _mm_load_si128((__m128i const *) ref_bases);

			// load relevant matrix cells (up, diag, left)
			__m128i up = _mm_load_si128((__m128i *) (scoreMat + SSE_SIZE * (prev_row_score * (refLength + 1) + ref_pos + 1)));
			__m128i diag = _mm_load_si128((__m128i *) (scoreMat + SSE_SIZE * (prev_row_score * (refLength + 1) + ref_pos)));
			__m128i left = _mm_load_si128((__m128i *) (scoreMat + SSE_SIZE * (current_row_score * (refLength + 1) + ref_pos)));

			// add gap penalties to up and left
			up = _mm_add_epi16(up, x_scoreGapRef);
			left = _mm_add_epi16(left, x_scoreGapRead);

			// UC ref
			sse_ref_bases = _mm_and_si128(sse_ref_bases,x_UCMask);

			__m128i valid_ref_base = _mm_or_si128(_mm_cmpeq_epi16(sse_ref_bases,x_C),_mm_or_si128(_mm_cmpeq_epi16(sse_ref_bases,x_G),_mm_or_si128(_mm_cmpeq_epi16(sse_ref_bases,x_T),_mm_cmpeq_epi16(sse_ref_bases,x_A))));

			//std::cout << "Ref base current:\t" << __m128i_toString<char>(sse_ref_bases);
			//std::cout << std::endl << "Valid:\t" << __m128i_toString<short>(valid_ref_base) << std::endl;

			__m128i valid_comp = _mm_and_si128(valid_read_base, valid_ref_base);

			// match read and ref bases
			// matches will hold 1, mismatches will hold 0
			__m128i match = _mm_cmpeq_epi16(sse_read_bases, sse_ref_bases);

			// bitwise and between 0 (mismatches) and match score gives match score only for equal bases
			// add these to diagonal value

			diag = _mm_add_epi16(diag, _mm_and_si128(valid_comp,_mm_and_si128(match, x_scoreMatch)));
			//diag = _mm_add_epi16(diag, _mm_and_si128(match, x_scoreMatch));
			// bitwise not and between 1 (matches) and mismatch score gives mismatch score only for unequal bases
			// add these to diagonal value
			diag = _mm_add_epi16(diag, _mm_and_si128(valid_comp,_mm_andnot_si128(match, x_scoreMismatch)));

			// Cell value will be max of upper + gap penalty, left + gap penalty, diag + match/mismatch score or 0
			__m128i cell = _mm_max_epi16(diag, _mm_max_epi16(left, _mm_max_epi16(up, x_zeros)));

			// Determine pointer direction

			__m128i pointer = x_zeros;

			// up
			match = _mm_cmpeq_epi16(cell, up);
			pointer = _mm_max_epi16(pointer, _mm_and_si128(match, x_p_up));

			// left
			match = _mm_cmpeq_epi16(cell, left);
			pointer = _mm_max_epi16(pointer, _mm_and_si128(match, x_p_left));

			// diag
			// Set diag pointer only if comparison of ATGC
			match =  _mm_and_si128(valid_comp,_mm_cmpeq_epi16(cell, diag));
			pointer = _mm_max_epi16(pointer, _mm_and_si128(match, x_p_diag));

			// Store score and pointer
			_mm_store_si128((__m128i *) (scoreMat + SSE_SIZE * (current_row_score * (refLength + 1) + ref_pos + 1)), cell);

			std::cout << __m128i_toString<short>(cell) << std::endl;

			_mm_store_si128((__m128i *) (matrix + SSE_SIZE * (current_row_aln * (refLength + 1) + ref_pos + 1)), pointer);

			// Store read and ref positions if new max score
			match = _mm_cmpgt_epi16(cell, max_score);

			best_read_pos = _mm_max_epi16(best_read_pos,_mm_and_si128(match, sse_read_pos));
			best_ref_pos = _mm_max_epi16(best_ref_pos,_mm_and_si128(match, sse_ref_pos));

			//std::cout << "Best read pos\t" << __m128i_toString<short>(best_read_pos);
			//std::cout << "Best ref pos\t" << __m128i_toString<short>(best_ref_pos);

			//char a;
			//std::cin >> a;

			max_score = _mm_max_epi16(cell, max_score);

			sse_ref_pos = _mm_add_epi16(sse_ref_pos, see_increment);

			//std::cout << "ref pos " << __m128i_toString<short>(sse_ref_pos);

		}

		std::cout << std::endl;
		prev_row_score = current_row_score;
		(++current_row_score) &= 1;

		++current_row_aln;

		sse_read_pos = _mm_add_epi16(sse_read_pos, see_increment);

		//std::cout << __m128i_toString<short>(sse_read_pos);

	}

	std::cout << "Max scores: " << __m128i_toString<short>(max_score);

	free(scoreMat);

	_mm_store_si128((__m128i *)best_coordinates, best_read_pos);
	_mm_store_si128((__m128i *)best_coordinates + 1, best_ref_pos);

}

void SSEKernel::calc_alignment(char const * const * const read,
				char const * const * const ref, Alignment * const alignment) {

	std::cout << "Aln Length:\t" << alnLength << std::endl;

	short * matrix = 0;
	malloc16(matrix, sizeof(short) * (refLength + 1) * (readLength + 1) * SSE_SIZE,16);
	memset(matrix, 0, (refLength + 1) * (readLength + 1) * SSE_SIZE * sizeof(short));

	short * best_coordinates = 0;
	malloc16(best_coordinates, sizeof(short) * 2 * SSE_SIZE,16);
	memset(best_coordinates, 0, sizeof(short) * 2 * SSE_SIZE);

	std::cout << "Score matrix" << std::endl;

	calc_alignment_matrix(read, ref, matrix, best_coordinates);

	for (int SSE_register = 0; SSE_register < SSE_SIZE; ++SSE_register) {
		std::cout << "Matrix:" << std::endl;
		for (int i = 0; i < readLength + 1; ++i) {
			for (int j = 0; j < refLength + 1; ++j) {
				//__m128i cell = _mm_load_si128((__m128i *) (matrix + SSE_SIZE * (i * (refLength + 1) + j)));
				//std::cout << __m128i_toString<short>(cell) << std::endl;
				std::cout << *(matrix + SSE_SIZE * (i * (refLength + 1) + j) + SSE_register) << " ";
			}
			std::cout << std::endl;
		}
	}

	char * alignments = new char[alnLength * 2 * SSE_SIZE];

	// Retreive best read and ref pos
	__m128i best_read_positions = _mm_load_si128((__m128i *) best_coordinates);
	__m128i best_ref_positions = _mm_load_si128((__m128i *) best_coordinates + 1);

	std::cout << "Best read: " << __m128i_toString<short>(best_read_positions);
	std::cout << std::endl << "Best ref: " << __m128i_toString<short>(best_ref_positions);
	std::cout << std::endl;

	for (int SSE_register = 0; SSE_register < SSE_SIZE; ++SSE_register) {
		std::cout << "register " << SSE_register << std::endl;
		short read_pos = *(best_coordinates + SSE_register);
		short ref_pos = *(best_coordinates + SSE_SIZE + SSE_register);
		std::cout << "Cur read_pos " << read_pos << " Cur ref_pos " << ref_pos << std::endl;

		int aln_pos = alnLength - 2;

		short backtrack = *(matrix + SSE_SIZE * ((read_pos + 1) * (refLength + 1) + ref_pos + 1) + SSE_register);

		while (backtrack != START) {
			std::cout << "Cur backtrack " << backtrack << std::endl;

			if (backtrack == UP) {
				alignments[(SSE_register * alnLength * 2) + alnLength + aln_pos] = '-';
				alignments[(SSE_register * alnLength * 2) + aln_pos] = read[SSE_register][read_pos--];
			}

			char base = *(read[SSE_register] + read_pos);
			std::cout << "Read base:\t" << base << std::endl;
			if (backtrack == LEFT) {
				alignments[(SSE_register * alnLength * 2) + aln_pos] = '-';
				alignments[(SSE_register * alnLength * 2) + alnLength + aln_pos] = ref[SSE_register][ref_pos--];
			}
			base = *(ref[SSE_register] + ref_pos);
			std::cout << "Ref base:\t" << base << std::endl;
			if (backtrack == DIAG) {
				alignments[(SSE_register * alnLength * 2) + aln_pos] = read[SSE_register][read_pos--];
				alignments[(SSE_register * alnLength * 2) + alnLength + aln_pos] = ref[SSE_register][ref_pos--];
			}

			backtrack = *(matrix + SSE_SIZE * ((read_pos + 1) * (refLength + 1) + ref_pos + 1) + SSE_register);
			--aln_pos;

			std::cout << "Cur read pos:\t" << read_pos << std::endl << "Cur refpos:\t" << ref_pos << std::endl;

		}

		alignment[SSE_register].read = new char[alnLength];
		alignment[SSE_register].ref = new char[alnLength];

		memcpy(alignment[SSE_register].read, alignments + (SSE_register * alnLength * 2), alnLength * sizeof(char));
		memcpy(alignment[SSE_register].ref, alignments + (SSE_register * alnLength * 2) + alnLength, alnLength * sizeof(char));

		alignment[SSE_register].readStart = aln_pos + 1;
		alignment[SSE_register].refStart = aln_pos + 1;

		alignment[SSE_register].readEnd = alnLength - 1;
		alignment[SSE_register].refEnd = alnLength - 1;

		std::cout << "==================" << std::endl << "\"";
		std::cout << alignment[SSE_register].read + alignment[SSE_register].readStart;
		std::cout << "\"" << std::endl << "\"";
		std::cout << alignment[SSE_register].ref + alignment[SSE_register].refStart;
		std::cout << "\"" << std::endl << "==================" << std::endl;
	}



//		alignment->ref = new char[alnLength];
//
//		memcpy(alignment->read, alignments, alnLength * sizeof(char));
//		memcpy(alignment->ref, alignments + alnLength, alnLength * sizeof(char));
//
//		alignment->readStart = aln_pos + 1;
//		alignment->refStart = aln_pos + 1;
//
//		alignment->readEnd = alnLength - 1;
//		alignment->refEnd = alnLength - 1;
//
//		std::cout << "\tBacktrack end" << std::endl;
}


void SSEKernel::score_alignment (char const * const * const read, char const * const * const ref, short * const scores) {

	__m128i max_score = x_zeros;

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

		// UC read
		sse_read_bases = _mm_and_si128(sse_read_bases,x_UCMask);

		//std::cout << "Read base:\t" << __m128i_toString<char>(sse_read_bases);

		__m128i valid_read_base = _mm_or_si128(_mm_cmpeq_epi16(sse_read_bases,x_C),_mm_or_si128(_mm_cmpeq_epi16(sse_read_bases,x_G),_mm_or_si128(_mm_cmpeq_epi16(sse_read_bases,x_T),_mm_cmpeq_epi16(sse_read_bases,x_A))));

		for (int ref_pos = 0; ref_pos < refLength; ++ref_pos) {

			// load ref base
			for (int base = 0; base < SSE_SIZE; ++base) {
				ref_bases[base] = *(ref[base] + ref_pos);
			}

			sse_ref_bases = _mm_load_si128((__m128i const *) ref_bases);

			// load relevant matrix cells (up, diag, left)
			__m128i up = _mm_load_si128((__m128i *) (matrix + SSE_SIZE * (prev_row * (refLength + 1) + ref_pos + 1)));
			__m128i diag = _mm_load_si128((__m128i *) (matrix + SSE_SIZE * (prev_row * (refLength + 1) + ref_pos)));
			__m128i left = _mm_load_si128((__m128i *) (matrix + SSE_SIZE * (cur_row * (refLength + 1) + ref_pos)));

			// add gap penalties to up and left
			up = _mm_add_epi16(up, x_scoreGapRef);
			left = _mm_add_epi16(left, x_scoreGapRead);

			// UC ref
			sse_ref_bases = _mm_and_si128(sse_ref_bases,x_UCMask);

			__m128i valid_ref_base = _mm_or_si128(_mm_cmpeq_epi16(sse_ref_bases,x_C),_mm_or_si128(_mm_cmpeq_epi16(sse_ref_bases,x_G),_mm_or_si128(_mm_cmpeq_epi16(sse_ref_bases,x_T),_mm_cmpeq_epi16(sse_ref_bases,x_A))));

			__m128i valid_comp = _mm_and_si128(valid_read_base, valid_ref_base);

			//std::cout << "Ref base:\t" << __m128i_toString<char>(sse_ref_bases);

			// match read and ref bases
			// matches will hold 1, mismatches will hold 0
			__m128i match = _mm_cmpeq_epi16(sse_read_bases, sse_ref_bases);

			// bitwise and between 0 (mismatches) and match score gives match score only for equal bases
			// add these to diagonal value
			diag = _mm_add_epi16(diag, _mm_and_si128(valid_comp,_mm_and_si128(match, x_scoreMatch)));
			//diag = _mm_add_epi16(diag, _mm_and_si128(match, x_scoreMatch));
			// bitwise not and between 1 (matches) and mismatch score gives mismatch score only for unequal bases
			// add these to diagonal value
			//diag = _mm_add_epi16(diag, _mm_andnot_si128(match, x_scoreMismatch));
			diag = _mm_add_epi16(diag, _mm_and_si128(valid_comp,_mm_andnot_si128(match, x_scoreMismatch)));

			// Cell value will be max of upper + gap penalty, left + gap penalty, diag + match/mismatch score or 0
			__m128i cell = _mm_max_epi16(diag, _mm_max_epi16(left, _mm_max_epi16(up, x_zeros)));

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
