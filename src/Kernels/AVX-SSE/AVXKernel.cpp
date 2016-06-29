/*
 * AVXKernel.cpp
 *
 *  Created on: May 23, 2016
 *      Author: tobias.neumann
 */

#include "AVXKernel.h"
#include <memory.h>
#include <iostream>
#include <cstdio>
#include <sstream>
#include "../../util/versalignUtil.h"

#define SCORING_ROWS 2


template <typename T>
std::string __m256i_toString(const __m256i var) {
    std::stringstream sstr;
    const T* values = (const T*) &var;
    if (sizeof(T) == 1) {
        for (unsigned int i = 0; i < sizeof(__m256i); i++) {
            sstr << (int) values[i] << " ";
        }
    } else {
        for (unsigned int i = 0; i < sizeof(__m256i) / sizeof(T); i++) {
            sstr << values[i] << " ";
        }
    }
    return sstr.str();
}

void AVXKernel::calc_alignment_matrix(char const * const * const read,
		char const * const * const ref, short * const matrix, short * const best_coordinates) {

	// Tracking best read and ref positions
	__m256i best_read_pos = x_zeros;
	__m256i best_ref_pos = x_zeros;
	__m256i max_score = x_zeros;

	// Scoring matrix indices
	short prev_row_score = 0;
	short current_row_score = 1;

	// Alignment matrix index
	short current_row_aln = 1;

	short * scoreMat = 0;

	malloc16(scoreMat, sizeof(short) * (refLength + 1) * SCORING_ROWS * AVX_SIZE,16);
	memset(scoreMat, 0, sizeof(short) * (refLength + 1) * SCORING_ROWS * AVX_SIZE);

	// SSE conversion array for current read and ref bases
	align16 short read_bases [AVX_SIZE];
	align16 short ref_bases  [AVX_SIZE];

	// holds current read and ref base for SSE instruction
	__m256i sse_read_bases;
	__m256i sse_ref_bases;

	// holds current read and ref positions
	__m256i sse_read_pos = x_zeros;
	__m256i sse_ref_pos = x_zeros;

	// holds 1 for increment
	__m256i see_increment = short_to_avx(1);

	for (int read_pos = 0; read_pos < readLength; ++read_pos) {

		sse_ref_pos = x_zeros;

		// load read base
		for (int base = 0; base < AVX_SIZE; ++base) {
			read_bases[base] = *(read[base] + read_pos);
		}

		sse_read_bases = _mm256_load_si256((__m256i const *) read_bases);

		// UC read
		sse_read_bases = _mm256_and_si256(sse_read_bases,x_UCMask);

		//std::cout << "Read base current:\t" << __m256i_toString<char>(sse_read_bases);

		__m256i valid_read_base = _mm256_or_si256(_mm256_cmpeq_epi16(sse_read_bases,x_C),_mm256_or_si256(_mm256_cmpeq_epi16(sse_read_bases,x_G),_mm256_or_si256(_mm256_cmpeq_epi16(sse_read_bases,x_T),_mm256_cmpeq_epi16(sse_read_bases,x_A))));

		//std::cout << std::endl << "Valid:\t" << __m256i_toString<short>(valid_read_base) << std::endl;

		for (int ref_pos = 0; ref_pos < refLength; ++ref_pos) {

			// load ref base
			for (int base = 0; base < AVX_SIZE; ++base) {
				ref_bases[base] = *(ref[base] + ref_pos);
			}

			sse_ref_bases = _mm256_load_si256((__m256i const *) ref_bases);

			// load relevant matrix cells (up, diag, left)
			__m256i up = _mm256_load_si256((__m256i *) (scoreMat + AVX_SIZE * (prev_row_score * (refLength + 1) + ref_pos + 1)));
			__m256i diag = _mm256_load_si256((__m256i *) (scoreMat + AVX_SIZE * (prev_row_score * (refLength + 1) + ref_pos)));
			__m256i left = _mm256_load_si256((__m256i *) (scoreMat + AVX_SIZE * (current_row_score * (refLength + 1) + ref_pos)));

			// add gap penalties to up and left
			up = _mm256_add_epi16(up, x_scoreGapRef);
			left = _mm256_add_epi16(left, x_scoreGapRead);

			// UC ref
			sse_ref_bases = _mm256_and_si256(sse_ref_bases,x_UCMask);

			__m256i valid_ref_base = _mm256_or_si256(_mm256_cmpeq_epi16(sse_ref_bases,x_C),_mm256_or_si256(_mm256_cmpeq_epi16(sse_ref_bases,x_G),_mm256_or_si256(_mm256_cmpeq_epi16(sse_ref_bases,x_T),_mm256_cmpeq_epi16(sse_ref_bases,x_A))));

			//std::cout << "Ref base current:\t" << __m256i_toString<char>(sse_ref_bases);
			//std::cout << std::endl << "Valid:\t" << __m256i_toString<short>(valid_ref_base) << std::endl;

			__m256i valid_comp = _mm256_and_si256(valid_read_base, valid_ref_base);

			// match read and ref bases
			// matches will hold 1, mismatches will hold 0
			__m256i match = _mm256_cmpeq_epi16(sse_read_bases, sse_ref_bases);

			// bitwise and between 0 (mismatches) and match score gives match score only for equal bases
			// add these to diagonal value

			diag = _mm256_add_epi16(diag, _mm256_and_si256(valid_comp,_mm256_and_si256(match, x_scoreMatch)));
			//diag = _mm256_add_epi16(diag, _mm256_and_si256(match, x_scoreMatch));
			// bitwise not and between 1 (matches) and mismatch score gives mismatch score only for unequal bases
			// add these to diagonal value
			diag = _mm256_add_epi16(diag, _mm256_and_si256(valid_comp,_mm256_andnot_si256(match, x_scoreMismatch)));

			// Cell value will be max of upper + gap penalty, left + gap penalty, diag + match/mismatch score or 0
			__m256i cell = _mm256_max_epi16(diag, _mm256_max_epi16(left, _mm256_max_epi16(up, x_zeros)));

			// Determine pointer direction

			__m256i pointer = x_zeros;

			// up
			match = _mm256_cmpeq_epi16(cell, up);
			pointer = _mm256_max_epi16(pointer, _mm256_and_si256(match, x_p_up));

			// left
			match = _mm256_cmpeq_epi16(cell, left);
			pointer = _mm256_max_epi16(pointer, _mm256_and_si256(match, x_p_left));

			// diag
			// Set diag pointer only if comparison of ATGC
			match =  _mm256_and_si256(valid_comp,_mm256_cmpeq_epi16(cell, diag));
			pointer = _mm256_max_epi16(pointer, _mm256_and_si256(match, x_p_diag));

			// Store score and pointer
			_mm256_store_si256((__m256i *) (scoreMat + AVX_SIZE * (current_row_score * (refLength + 1) + ref_pos + 1)), cell);
			_mm256_store_si256((__m256i *) (matrix + AVX_SIZE * (current_row_aln * (refLength + 1) + ref_pos + 1)), pointer);

			// Store read and ref positions if new max score
			match = _mm256_cmpgt_epi16(cell, max_score);

//			std::cout << "Scores:\t" << __m256i_toString<short>(cell) << std::endl;
//			std::cout << "Max:\t" << __m256i_toString<short>(max_score) << std::endl;
//			std::cout << "Match:\t" << __m256i_toString<short>(match) << std::endl;

			best_read_pos = _mm256_max_epi16(_mm256_andnot_si256(match, best_read_pos),_mm256_and_si256(match, sse_read_pos));
			best_ref_pos = _mm256_max_epi16(_mm256_andnot_si256(match, best_ref_pos),_mm256_and_si256(match, sse_ref_pos));

			//std::cout << "Best read pos\t" << __m256i_toString<short>(best_read_pos);
			//std::cout << "Best ref pos\t" << __m256i_toString<short>(best_ref_pos);

			max_score = _mm256_max_epi16(cell, max_score);

			sse_ref_pos = _mm256_add_epi16(sse_ref_pos, see_increment);

			//std::cout << "ref pos " << __m256i_toString<short>(sse_ref_pos);

		}
		prev_row_score = current_row_score;
		(++current_row_score) &= 1;

		++current_row_aln;

		sse_read_pos = _mm256_add_epi16(sse_read_pos, see_increment);

		//std::cout << __m256i_toString<short>(sse_read_pos);

	}

//	std::cout << "Max scores: " << __m256i_toString<short>(max_score);

	free(scoreMat);

	_mm256_store_si256((__m256i *)best_coordinates, best_read_pos);
	_mm256_store_si256((__m256i *)best_coordinates + 1, best_ref_pos);

}

void AVXKernel::calculate_alignment_matrix_needleman_wunsch(char const * const * const read,
		char const * const * const ref, short * const matrix, short * const best_coordinates) {

	// Tracking positions where invalid characters start
	__m256i max_read_pos = short_to_avx(readLength - 1);
	__m256i max_ref_pos = short_to_avx(refLength - 1);

	// Tracking rowwise maxima
	__m256i rowMax = short_to_avx(SHRT_MIN);
	__m256i rowMaxIndex = x_zeros;

	// Tracking global maxima
//	__m256i globalRowMax = short_to_avx(SHRT_MIN);
	__m256i globalRowMaxIndex = short_to_avx(-1);

	// Scoring matrix indices
	short prev_row_score = 0;
	short current_row_score = 1;

	// Alignment matrix index
	short current_row_aln = 1;

	short * scoreMat = 0;

	malloc16(scoreMat, sizeof(short) * (refLength + 1) * SCORING_ROWS * AVX_SIZE,16);
	memset(scoreMat, 0, sizeof(short) * (refLength + 1) * SCORING_ROWS * AVX_SIZE);

	// SSE conversion array for current read and ref bases
	align16 short read_bases [AVX_SIZE];
	align16 short ref_bases  [AVX_SIZE];

	// holds current read and ref base for SSE instruction
	__m256i sse_read_bases;
	__m256i sse_ref_bases;

	// holds current read and ref positions
	__m256i sse_read_pos = x_zeros;
	__m256i sse_ref_pos = x_zeros;

	// holds 1 for increment
	__m256i see_increment = short_to_avx(1);

	for (int read_pos = 0; read_pos < readLength; ++read_pos) {

		sse_ref_pos = x_zeros;

		// load read base
		for (int base = 0; base < AVX_SIZE; ++base) {
			read_bases[base] = *(read[base] + read_pos);
		}

		sse_read_bases = _mm256_load_si256((__m256i const *) read_bases);

		// UC read
		sse_read_bases = _mm256_and_si256(sse_read_bases,x_UCMask);

		__m256i valid_read_base = _mm256_or_si256(_mm256_cmpeq_epi16(sse_read_bases,x_C),_mm256_or_si256(_mm256_cmpeq_epi16(sse_read_bases,x_G),_mm256_or_si256(_mm256_cmpeq_epi16(sse_read_bases,x_T),_mm256_cmpeq_epi16(sse_read_bases,x_A))));

		// First alignment column is always UP
		_mm256_store_si256((__m256i *) (matrix + AVX_SIZE * (current_row_aln * (refLength + 1))), x_p_up);
		// First score column is always continued gap-ref score
		_mm256_store_si256((__m256i *) (scoreMat + AVX_SIZE * (current_row_score * (refLength + 1))), short_to_avx((read_pos + 1) * scoreGapRef));

		// if max_readpos == readLength - 1 AND read char is invalid -> max_read_pos = read_pos - 1
		max_read_pos = _mm_blendv_si256(max_read_pos, _mm256_sub_epi16(sse_read_pos,see_increment), _mm256_andnot_si256(valid_read_base,_mm256_cmpeq_epi16(max_read_pos, short_to_avx(readLength - 1))));

//		std::cout << "Max read pos:\t" << __m256i_toString<short>(max_read_pos) << std::endl;
//		std::cout << "Read pos:\t" << __m256i_toString<short>(sse_read_pos) << std::endl;
//		std::cout << "Read Length:\t" << __m256i_toString<short>(short_to_avx(readLength - 1)) << std::endl;
//		std::cout << "Valid read base:\t" << __m256i_toString<short>(valid_read_base) << std::endl;
//		std::cout << "Max read pos == read Length - 1:\t" << __m256i_toString<short>(_mm256_cmpeq_epi16(max_read_pos, short_to_avx(readLength - 1))) << std::endl;
//		std::cout << "cmp:\t" << __m256i_toString<short>(_mm256_andnot_si256(valid_read_base,_mm256_cmpeq_epi16(max_read_pos, short_to_avx(readLength - 1)))) << std::endl;
//		std::cout << std::endl;

		// Save previous row max if read ends prematurely
		// if max_readpos + 1 == read_pos -> globalRowMax = rowMax; globalRowMaxIndex = rowMaxIndex;
//		__m256i sel = _mm256_cmpeq_epi16(_mm256_add_epi16(max_read_pos,see_increment), sse_read_pos);
//		globalRowMax =  _mm_blendv_si256(globalRowMax,rowMax,sel);
//		globalRowMaxIndex =  _mm_blendv_si256(globalRowMaxIndex,rowMaxIndex,sel);
		globalRowMaxIndex =  _mm_blendv_si256(globalRowMaxIndex,rowMaxIndex,_mm256_cmpeq_epi16(_mm256_add_epi16(max_read_pos,see_increment), sse_read_pos));

//		std::cout << "Max_readpos + 1 == readpos:\t" << __m256i_toString<short>(sel) << std::endl;
//		std::cout << "globalRowMax:\t" << __m256i_toString<short>(globalRowMax) << std::endl;
//		std::cout << "globalRowMaxIndex:\t" << __m256i_toString<short>(globalRowMaxIndex) << std::endl;
//		std::cout << "rowMax:\t" << __m256i_toString<short>(rowMax) << std::endl;
//		std::cout << "rowMaxIndex:\t" << __m256i_toString<short>(rowMaxIndex) << std::endl;

		rowMax = _mm256_load_si256((__m256i *) (scoreMat + AVX_SIZE * (current_row_score * (refLength + 1))));
		rowMaxIndex = x_zeros;

		for (int ref_pos = 0; ref_pos < refLength; ++ref_pos) {

			// load ref base
			for (int base = 0; base < AVX_SIZE; ++base) {
				ref_bases[base] = *(ref[base] + ref_pos);
			}

			sse_ref_bases = _mm256_load_si256((__m256i const *) ref_bases);

			// load relevant matrix cells (up, diag, left)
			__m256i up = _mm256_load_si256((__m256i *) (scoreMat + AVX_SIZE * (prev_row_score * (refLength + 1) + ref_pos + 1)));
			__m256i diag = _mm256_load_si256((__m256i *) (scoreMat + AVX_SIZE * (prev_row_score * (refLength + 1) + ref_pos)));
			__m256i left = _mm256_load_si256((__m256i *) (scoreMat + AVX_SIZE * (current_row_score * (refLength + 1) + ref_pos)));

			// add gap penalties to up and left
			up = _mm256_add_epi16(up, x_scoreGapRef);
			left = _mm256_add_epi16(left, x_scoreGapRead);

			// UC ref
			sse_ref_bases = _mm256_and_si256(sse_ref_bases,x_UCMask);

			__m256i valid_ref_base = _mm256_or_si256(_mm256_cmpeq_epi16(sse_ref_bases,x_C),_mm256_or_si256(_mm256_cmpeq_epi16(sse_ref_bases,x_G),_mm256_or_si256(_mm256_cmpeq_epi16(sse_ref_bases,x_T),_mm256_cmpeq_epi16(sse_ref_bases,x_A))));

//			std::cout << "Ref base current:\t";
//			print_sse_char(sse_ref_bases);
//			std::cout << std::endl << "Valid:\t" << __m256i_toString<short>(valid_ref_base) << std::endl;

			__m256i valid_comp = _mm256_and_si256(valid_read_base, valid_ref_base);

			// match read and ref bases
			// matches will hold 1, mismatches will hold 0
			__m256i match = _mm256_cmpeq_epi16(sse_read_bases, sse_ref_bases);

			// bitwise and between 0 (mismatches) and match score gives match score only for equal bases
			// add these to diagonal value

			diag = _mm256_add_epi16(diag, _mm256_and_si256(valid_comp,_mm256_and_si256(match, x_scoreMatch)));
			//diag = _mm256_add_epi16(diag, _mm256_and_si256(match, x_scoreMatch));
			// bitwise not and between 1 (matches) and mismatch score gives mismatch score only for unequal bases
			// add these to diagonal value
			diag = _mm256_add_epi16(diag, _mm256_and_si256(valid_comp,_mm256_andnot_si256(match, x_scoreMismatch)));

			// Cell value will be max of upper + gap penalty, left + gap penalty, diag + match/mismatch score or 0
			__m256i cell = _mm256_max_epi16(diag, _mm256_max_epi16(left, up));

			// Determine pointer direction

			__m256i pointer = x_zeros;

			// up
			match = _mm256_cmpeq_epi16(cell, up);
			pointer = _mm256_max_epi16(pointer, _mm256_and_si256(match, x_p_up));

			// left
			match = _mm256_cmpeq_epi16(cell, left);
			pointer = _mm256_max_epi16(pointer, _mm256_and_si256(match, x_p_left));

			// diag
			// Set diag pointer only if comparison of ATGC
			match =  _mm256_and_si256(valid_comp,_mm256_cmpeq_epi16(cell, diag));
			pointer = _mm256_max_epi16(pointer, _mm256_and_si256(match, x_p_diag));

			// Store score and pointer
			_mm256_store_si256((__m256i *) (scoreMat + AVX_SIZE * (current_row_score * (refLength + 1) + ref_pos + 1)), cell);
			_mm256_store_si256((__m256i *) (matrix + AVX_SIZE * (current_row_aln * (refLength + 1) + ref_pos + 1)), pointer);

			// Store read and ref positions if new max score
			//match = _mm256_cmpgt_epi16(cell, max_score);

//			std::cout << "Scores:\t" << __m256i_toString<short>(cell) << std::endl;
//			std::cout << "Max:\t" << __m256i_toString<short>(max_score) << std::endl;
//			std::cout << "Match:\t" << __m256i_toString<short>(match) << std::endl;

			//best_read_pos = _mm256_max_epi16(_mm256_andnot_si256(match, best_read_pos),_mm256_and_si256(match, sse_read_pos));
			//best_ref_pos = _mm256_max_epi16(_mm256_andnot_si256(match, best_ref_pos),_mm256_and_si256(match, sse_ref_pos));

			//std::cout << "Best read pos\t" << __m256i_toString<short>(best_read_pos);
			//std::cout << "Best ref pos\t" << __m256i_toString<short>(best_ref_pos);

			// if max_ref_pos == refLength - 1 AND ref char is invalid -> max_ref_pos = ref_pos - 1
//			std::cout << std::endl << "Max_ref_pos:\t" << __m256i_toString<short>(max_ref_pos) << std::endl;
//			std::cout << std::endl << "sse_ref_pos:\t" << __m256i_toString<short>(sse_ref_pos) << std::endl;
//			std::cout << std::endl << "reflength - 1:\t" << __m256i_toString<short>(short_to_avx(refLength - 1)) << std::endl;
//			std::cout << std::endl << "valid_ref_base:\t" << __m256i_toString<short>(valid_ref_base) << std::endl;
//
//			char a;
//			std::cin >> a;
			max_ref_pos = _mm_blendv_si256(max_ref_pos, _mm256_sub_epi16(sse_ref_pos,see_increment), _mm256_andnot_si256(valid_ref_base,_mm256_cmpeq_epi16(max_ref_pos, short_to_avx(refLength - 1))));

			// if cur > rowMax => rowMax = cur; rowMaxIndex = ref_pos;
			__m256i sel = _mm256_cmpgt_epi16(cell, rowMax);
			rowMax = _mm_blendv_si256(rowMax,cell,sel);
			rowMaxIndex = _mm_blendv_si256(rowMaxIndex,sse_ref_pos,sel);

			sse_ref_pos = _mm256_add_epi16(sse_ref_pos, see_increment);

			//std::cout << "ref pos " << __m256i_toString<short>(sse_ref_pos);

		}
		prev_row_score = current_row_score;
		(++current_row_score) &= 1;

		++current_row_aln;

		sse_read_pos = _mm256_add_epi16(sse_read_pos, see_increment);

		//std::cout << __m256i_toString<short>(sse_read_pos);

	}

//	std::cout << "Max scores: " << __m256i_toString<short>(max_score);

	free(scoreMat);

//	__m256i tmp = _mm256_cmpgt_epi16(x_zeros,globalRowMaxIndex);

//	globalRowMax = _mm_blendv_si256(globalRowMax, rowMax, tmp);
//	globalRowMaxIndex = _mm_blendv_si256(globalRowMaxIndex, rowMaxIndex, tmp);
	globalRowMaxIndex = _mm_blendv_si256(globalRowMaxIndex, rowMaxIndex, _mm256_cmpgt_epi16(x_zeros,globalRowMaxIndex));

//	std::cout << "Max_ref_pos:\t" << __m256i_toString<short>(max_ref_pos) << std::endl;
//	std::cout << "globalRowMaxIndex:\t" << __m256i_toString<short>(globalRowMaxIndex) << std::endl;

	__m256i best_ref_pos = _mm256_min_epi16(max_ref_pos,globalRowMaxIndex);

	_mm256_store_si256((__m256i *)best_coordinates, max_read_pos);
	_mm256_store_si256((__m256i *)best_coordinates + 1, best_ref_pos);
}

void AVXKernel::calc_alignment(char const * const * const read,
				char const * const * const ref, Alignment * const alignment) {

	std::cout << "Aln Length:\t" << alnLength << std::endl;

	short * matrix = 0;
	malloc16(matrix, sizeof(short) * (refLength + 1) * (readLength + 1) * AVX_SIZE,16);
	memset(matrix, 0, (refLength + 1) * (readLength + 1) * AVX_SIZE * sizeof(short));

	short * best_coordinates = 0;
	malloc16(best_coordinates, sizeof(short) * 2 * AVX_SIZE,16);
	memset(best_coordinates, 0, sizeof(short) * 2 * AVX_SIZE);

	std::cout << "Score matrix" << std::endl;

	calc_alignment_matrix(read, ref, matrix, best_coordinates);

	for (int SSE_register = 0; SSE_register < AVX_SIZE; ++SSE_register) {
		std::cout << "Matrix:" << std::endl;
		for (int i = 0; i < readLength + 1; ++i) {
			for (int j = 0; j < refLength + 1; ++j) {
				//__m256i cell = _mm256_load_si256((__m256i *) (matrix + AVX_SIZE * (i * (refLength + 1) + j)));
				//std::cout << __m256i_toString<short>(cell) << std::endl;
				std::cout << *(matrix + AVX_SIZE * (i * (refLength + 1) + j) + SSE_register) << " ";
			}
			std::cout << std::endl;
		}
	}

	char * alignments = new char[alnLength * 2 * AVX_SIZE];

	// Retreive best read and ref pos
	__m256i best_read_positions = _mm256_load_si256((__m256i *) best_coordinates);
	__m256i best_ref_positions = _mm256_load_si256((__m256i *) best_coordinates + 1);

//	std::cout << "Best read: " << __m256i_toString<short>(best_read_positions);
//	std::cout << std::endl << "Best ref: " << __m256i_toString<short>(best_ref_positions);
//	std::cout << std::endl;

	for (int SSE_register = 0; SSE_register < AVX_SIZE; ++SSE_register) {
//		std::cout << "register " << SSE_register << std::endl;
		short read_pos = *(best_coordinates + SSE_register);
		short ref_pos = *(best_coordinates + AVX_SIZE + SSE_register);
//		std::cout << "Cur read_pos " << read_pos << " Cur ref_pos " << ref_pos << std::endl;

		int aln_pos = alnLength - 2;

		alignments[(SSE_register * alnLength * 2) + alnLength - 1] = '\0';
		alignments[(SSE_register * alnLength * 2) + (2 * alnLength) - 1] = '\0';

		short backtrack = *(matrix + AVX_SIZE * ((read_pos + 1) * (refLength + 1) + ref_pos + 1) + SSE_register);

		while (backtrack != START) {
//			std::cout << "Cur backtrack " << backtrack << std::endl;

			if (backtrack == UP) {
				alignments[(SSE_register * alnLength * 2) + alnLength + aln_pos] = '-';
				alignments[(SSE_register * alnLength * 2) + aln_pos] = read[SSE_register][read_pos--];
			}

//			char base = *(read[SSE_register] + read_pos);
//			std::cout << "Read base:\t" << base << std::endl;
			if (backtrack == LEFT) {
				alignments[(SSE_register * alnLength * 2) + aln_pos] = '-';
				alignments[(SSE_register * alnLength * 2) + alnLength + aln_pos] = ref[SSE_register][ref_pos--];
			}
//			base = *(ref[SSE_register] + ref_pos);
//			std::cout << "Ref base:\t" << base << std::endl;
			if (backtrack == DIAG) {
				alignments[(SSE_register * alnLength * 2) + aln_pos] = read[SSE_register][read_pos--];
				alignments[(SSE_register * alnLength * 2) + alnLength + aln_pos] = ref[SSE_register][ref_pos--];
			}

			backtrack = *(matrix + AVX_SIZE * ((read_pos + 1) * (refLength + 1) + ref_pos + 1) + SSE_register);
			--aln_pos;

//			std::cout << "Cur read pos:\t" << read_pos << std::endl << "Cur refpos:\t" << ref_pos << std::endl;

		}

		alignment[SSE_register].read = new char[alnLength];
		alignment[SSE_register].ref = new char[alnLength];

		memcpy(alignment[SSE_register].read, alignments + (SSE_register * alnLength * 2), alnLength * sizeof(char));
		memcpy(alignment[SSE_register].ref, alignments + (SSE_register * alnLength * 2) + alnLength, alnLength * sizeof(char));

		alignment[SSE_register].readStart = aln_pos + 1;
		alignment[SSE_register].refStart = aln_pos + 1;

		alignment[SSE_register].readEnd = alnLength - 1;
		alignment[SSE_register].refEnd = alnLength - 1;
	}

	delete []alignments; alignments = 0;
	free(best_coordinates);
	free(matrix);
}

void AVXKernel::calc_alignment_needleman_wunsch(char const * const * const read,
		char const * const * const ref, Alignment * const alignment) {


	std::cout << "Aln Length:\t" << alnLength << std::endl;

	short * matrix = 0;
	malloc16(matrix, sizeof(short) * (refLength + 1) * (readLength + 1) * AVX_SIZE,16);
	memset(matrix, 0, (refLength + 1) * (readLength + 1) * AVX_SIZE * sizeof(short));

	short * best_coordinates = 0;
	malloc16(best_coordinates, sizeof(short) * 2 * AVX_SIZE,16);
	memset(best_coordinates, 0, sizeof(short) * 2 * AVX_SIZE);

	//std::cout << "Score matrix" << std::endl;

	calculate_alignment_matrix_needleman_wunsch(read, ref, matrix, best_coordinates);

	for (int SSE_register = 0; SSE_register < AVX_SIZE; ++SSE_register) {
		//std::cout << "Matrix:" << std::endl;
		for (int i = 0; i < readLength + 1; ++i) {
			for (int j = 0; j < refLength + 1; ++j) {
				//__m256i cell = _mm256_load_si256((__m256i *) (matrix + AVX_SIZE * (i * (refLength + 1) + j)));
				//std::cout << __m256i_toString<short>(cell) << std::endl;
				//std::cout << *(matrix + AVX_SIZE * (i * (refLength + 1) + j) + SSE_register) << " ";
			}
			//std::cout << std::endl;
		}
	}

	char * alignments = new char[alnLength * 2 * AVX_SIZE];

	// Retreive best read and ref pos
	__m256i best_read_positions = _mm256_load_si256((__m256i *) best_coordinates);
	__m256i best_ref_positions = _mm256_load_si256((__m256i *) best_coordinates + 1);

	std::cout << "Best read: " << __m256i_toString<short>(best_read_positions);
	std::cout << std::endl << "Best ref: " << __m256i_toString<short>(best_ref_positions);
	std::cout << std::endl;

	for (int SSE_register = 0; SSE_register < AVX_SIZE; ++SSE_register) {
		//		std::cout << "register " << SSE_register << std::endl;
		short read_pos = *(best_coordinates + SSE_register);
		short ref_pos = *(best_coordinates + AVX_SIZE + SSE_register);
		//		std::cout << "Cur read_pos " << read_pos << " Cur ref_pos " << ref_pos << std::endl;

		int aln_pos = alnLength - 2;

		alignments[(SSE_register * alnLength * 2) + alnLength - 1] = '\0';
		alignments[(SSE_register * alnLength * 2) + (2 * alnLength) - 1] = '\0';

		short backtrack = *(matrix + AVX_SIZE * ((read_pos + 1) * (refLength + 1) + ref_pos + 1) + SSE_register);

		while (backtrack != START) {
			//			std::cout << "Cur backtrack " << backtrack << std::endl;

			if (backtrack == UP) {
				alignments[(SSE_register * alnLength * 2) + alnLength + aln_pos] = '-';
				alignments[(SSE_register * alnLength * 2) + aln_pos] = read[SSE_register][read_pos--];
			}

			//			char base = *(read[SSE_register] + read_pos);
			//			std::cout << "Read base:\t" << base << std::endl;
			if (backtrack == LEFT) {
				alignments[(SSE_register * alnLength * 2) + aln_pos] = '-';
				alignments[(SSE_register * alnLength * 2) + alnLength + aln_pos] = ref[SSE_register][ref_pos--];
			}
			//			base = *(ref[SSE_register] + ref_pos);
			//			std::cout << "Ref base:\t" << base << std::endl;
			if (backtrack == DIAG) {
				alignments[(SSE_register * alnLength * 2) + aln_pos] = read[SSE_register][read_pos--];
				alignments[(SSE_register * alnLength * 2) + alnLength + aln_pos] = ref[SSE_register][ref_pos--];
			}

			backtrack = *(matrix + AVX_SIZE * ((read_pos + 1) * (refLength + 1) + ref_pos + 1) + SSE_register);
			--aln_pos;

			//			std::cout << "Cur read pos:\t" << read_pos << std::endl << "Cur refpos:\t" << ref_pos << std::endl;

		}

		alignment[SSE_register].read = new char[alnLength];
		alignment[SSE_register].ref = new char[alnLength];

		memcpy(alignment[SSE_register].read, alignments + (SSE_register * alnLength * 2), alnLength * sizeof(char));
		memcpy(alignment[SSE_register].ref, alignments + (SSE_register * alnLength * 2) + alnLength, alnLength * sizeof(char));

		alignment[SSE_register].readStart = aln_pos + 1;
		alignment[SSE_register].refStart = aln_pos + 1;

		alignment[SSE_register].readEnd = alnLength - 1;
		alignment[SSE_register].refEnd = alnLength - 1;
	}

	delete []alignments; alignments = 0;
	free(best_coordinates);
	free(matrix);

}


void AVXKernel::score_alignment (char const * const * const read, char const * const * const ref, short * const scores) {

	__m256i max_score = x_zeros;

	// Initialize SSE matrix
	short * matrix = 0;

	malloc16(matrix, sizeof(short) * (refLength + 1) * SCORING_ROWS * AVX_SIZE,16);
	memset(matrix, 0, sizeof(short) * (refLength + 1) * SCORING_ROWS * AVX_SIZE);

	// offset for first and second row
	int prev_row = 0;
	int cur_row = 1;

	// SSE conversion array for current read and ref bases
	align16 short read_bases [AVX_SIZE];
	align16 short ref_bases  [AVX_SIZE];

	// holds current read and ref base for SSE instruction
	__m256i sse_read_bases;
	__m256i sse_ref_bases;

	for (int read_pos = 0; read_pos < readLength; ++read_pos) {

		// load read base
		for (int base = 0; base < AVX_SIZE; ++base) {
			read_bases[base] = *(read[base] + read_pos);
		}

		sse_read_bases = _mm256_load_si256((__m256i const *) read_bases);

		// UC read
		sse_read_bases = _mm256_and_si256(sse_read_bases,x_UCMask);

		//std::cout << "Read base:\t" << __m256i_toString<char>(sse_read_bases);

		__m256i valid_read_base = _mm256_or_si256(_mm256_cmpeq_epi16(sse_read_bases,x_C),_mm256_or_si256(_mm256_cmpeq_epi16(sse_read_bases,x_G),_mm256_or_si256(_mm256_cmpeq_epi16(sse_read_bases,x_T),_mm256_cmpeq_epi16(sse_read_bases,x_A))));

		for (int ref_pos = 0; ref_pos < refLength; ++ref_pos) {

			// load ref base
			for (int base = 0; base < AVX_SIZE; ++base) {
				ref_bases[base] = *(ref[base] + ref_pos);
			}

			sse_ref_bases = _mm256_load_si256((__m256i const *) ref_bases);

			// load relevant matrix cells (up, diag, left)
			__m256i up = _mm256_load_si256((__m256i *) (matrix + AVX_SIZE * (prev_row * (refLength + 1) + ref_pos + 1)));
			__m256i diag = _mm256_load_si256((__m256i *) (matrix + AVX_SIZE * (prev_row * (refLength + 1) + ref_pos)));
			__m256i left = _mm256_load_si256((__m256i *) (matrix + AVX_SIZE * (cur_row * (refLength + 1) + ref_pos)));

			// add gap penalties to up and left
			up = _mm256_add_epi16(up, x_scoreGapRef);
			left = _mm256_add_epi16(left, x_scoreGapRead);

			// UC ref
			sse_ref_bases = _mm256_and_si256(sse_ref_bases,x_UCMask);

			__m256i valid_ref_base = _mm256_or_si256(_mm256_cmpeq_epi16(sse_ref_bases,x_C),_mm256_or_si256(_mm256_cmpeq_epi16(sse_ref_bases,x_G),_mm256_or_si256(_mm256_cmpeq_epi16(sse_ref_bases,x_T),_mm256_cmpeq_epi16(sse_ref_bases,x_A))));

			__m256i valid_comp = _mm256_and_si256(valid_read_base, valid_ref_base);

			//std::cout << "Ref base:\t" << __m256i_toString<char>(sse_ref_bases);

			// match read and ref bases
			// matches will hold 1, mismatches will hold 0
			__m256i match = _mm256_cmpeq_epi16(sse_read_bases, sse_ref_bases);

			// bitwise and between 0 (mismatches) and match score gives match score only for equal bases
			// add these to diagonal value
			diag = _mm256_add_epi16(diag, _mm256_and_si256(valid_comp,_mm256_and_si256(match, x_scoreMatch)));
			//diag = _mm256_add_epi16(diag, _mm256_and_si256(match, x_scoreMatch));
			// bitwise not and between 1 (matches) and mismatch score gives mismatch score only for unequal bases
			// add these to diagonal value
			//diag = _mm256_add_epi16(diag, _mm256_andnot_si256(match, x_scoreMismatch));
			diag = _mm256_add_epi16(diag, _mm256_and_si256(valid_comp,_mm256_andnot_si256(match, x_scoreMismatch)));

			// Cell value will be max of upper + gap penalty, left + gap penalty, diag + match/mismatch score or 0
			__m256i cell = _mm256_max_epi16(diag, _mm256_max_epi16(left, _mm256_max_epi16(up, x_zeros)));

			_mm256_store_si256((__m256i *) (matrix + AVX_SIZE * (cur_row * (refLength + 1) + ref_pos + 1)), cell);

			max_score = _mm256_max_epi16(cell, max_score);

		}
		prev_row = cur_row;
		(++cur_row) &= 1;
	}

	free(matrix);

	_mm256_store_si256((__m256i *)scores, max_score);
}


void AVXKernel::score_alignment_needleman_wunsch(char const * const * const read,
		char const * const * const ref, short * const scores) {

	__m256i max_score = x_zeros;

	// Initialize SSE matrix
	short * matrix = 0;

	malloc16(matrix, sizeof(short) * (refLength + 1) * SCORING_ROWS * AVX_SIZE,16);
	memset(matrix, 0, sizeof(short) * (refLength + 1) * SCORING_ROWS * AVX_SIZE);

	// offset for first and second row
	int prev_row = 0;
	int cur_row = 1;

	// SSE conversion array for current read and ref bases
	align16 short read_bases [AVX_SIZE];
	align16 short ref_bases  [AVX_SIZE];

	// holds current read and ref base for SSE instruction
	__m256i sse_read_bases;
	__m256i sse_ref_bases;

	for (int read_pos = 0; read_pos < readLength; ++read_pos) {

		// load read base
		for (int base = 0; base < AVX_SIZE; ++base) {
			read_bases[base] = *(read[base] + read_pos);
		}

		sse_read_bases = _mm256_load_si256((__m256i const *) read_bases);

		// UC read
		sse_read_bases = _mm256_and_si256(sse_read_bases,x_UCMask);

		//std::cout << "Read base:\t" << __m256i_toString<char>(sse_read_bases);

		__m256i valid_read_base = _mm256_or_si256(_mm256_cmpeq_epi16(sse_read_bases,x_C),_mm256_or_si256(_mm256_cmpeq_epi16(sse_read_bases,x_G),_mm256_or_si256(_mm256_cmpeq_epi16(sse_read_bases,x_T),_mm256_cmpeq_epi16(sse_read_bases,x_A))));

		for (int ref_pos = 0; ref_pos < refLength; ++ref_pos) {

			// load ref base
			for (int base = 0; base < AVX_SIZE; ++base) {
				ref_bases[base] = *(ref[base] + ref_pos);
			}

			sse_ref_bases = _mm256_load_si256((__m256i const *) ref_bases);

			// load relevant matrix cells (up, diag, left)
			__m256i up = _mm256_load_si256((__m256i *) (matrix + AVX_SIZE * (prev_row * (refLength + 1) + ref_pos + 1)));
			__m256i diag = _mm256_load_si256((__m256i *) (matrix + AVX_SIZE * (prev_row * (refLength + 1) + ref_pos)));
			__m256i left = _mm256_load_si256((__m256i *) (matrix + AVX_SIZE * (cur_row * (refLength + 1) + ref_pos)));

			// add gap penalties to up and left
			up = _mm256_add_epi16(up, x_scoreGapRef);
			left = _mm256_add_epi16(left, x_scoreGapRead);

			// UC ref
			sse_ref_bases = _mm256_and_si256(sse_ref_bases,x_UCMask);

			__m256i valid_ref_base = _mm256_or_si256(_mm256_cmpeq_epi16(sse_ref_bases,x_C),_mm256_or_si256(_mm256_cmpeq_epi16(sse_ref_bases,x_G),_mm256_or_si256(_mm256_cmpeq_epi16(sse_ref_bases,x_T),_mm256_cmpeq_epi16(sse_ref_bases,x_A))));

			__m256i valid_comp = _mm256_and_si256(valid_read_base, valid_ref_base);

			//std::cout << "Ref base:\t" << __m256i_toString<char>(sse_ref_bases);

			// match read and ref bases
			// matches will hold 1, mismatches will hold 0
			__m256i match = _mm256_cmpeq_epi16(sse_read_bases, sse_ref_bases);

			// bitwise and between 0 (mismatches) and match score gives match score only for equal bases
			// add these to diagonal value
			diag = _mm256_add_epi16(diag, _mm256_and_si256(valid_comp,_mm256_and_si256(match, x_scoreMatch)));
			//diag = _mm256_add_epi16(diag, _mm256_and_si256(match, x_scoreMatch));
			// bitwise not and between 1 (matches) and mismatch score gives mismatch score only for unequal bases
			// add these to diagonal value
			//diag = _mm256_add_epi16(diag, _mm256_andnot_si256(match, x_scoreMismatch));
			diag = _mm256_add_epi16(diag, _mm256_and_si256(valid_comp,_mm256_andnot_si256(match, x_scoreMismatch)));

			// Cell value will be max of upper + gap penalty, left + gap penalty, diag + match/mismatch score or 0
			__m256i cell = _mm256_max_epi16(diag, _mm256_max_epi16(left, up));

			_mm256_store_si256((__m256i *) (matrix + AVX_SIZE * (cur_row * (refLength + 1) + ref_pos + 1)), cell);
		}

		max_score = _mm256_max_epi16(max_score, _mm256_load_si256((__m256i *) (matrix + AVX_SIZE * (cur_row * (refLength + 1) + refLength))));

		prev_row = cur_row;
		(++cur_row) &= 1;
	}

	for (int ref_pos = 0; ref_pos < refLength + 1; ++ref_pos) {
		max_score = _mm256_max_epi16(max_score, _mm256_load_si256((__m256i *) (matrix + AVX_SIZE * (prev_row * (refLength + 1) + ref_pos))));
	}

	free(matrix);

	_mm256_store_si256((__m256i *)scores, max_score);
}

#undef SCORING_ROWS
#undef AVX_SIZE
