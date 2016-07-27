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
#include <string>

#define SCORING_ROWS 2

using std::string;
using std::to_string;

#ifndef NDEBUG

#include <sstream>

template<typename T>
std::string __m128i_toString(const __m128i var) {
	std::stringstream sstr;
	const T* values = (const T*) &var;
	if (sizeof(T) == 1) {
		for (unsigned int i = 0; i < sizeof(__m128i ); i++) {
			sstr << (int) values[i] << " ";
		}
	} else {
		for (unsigned int i = 0; i < sizeof(__m128i) / sizeof(T); i++) {
			sstr << values[i] << " ";
		}
	}
	return sstr.str();
}

#endif

void SSEKernel::compute_alignments(int const & opt, int const & aln_number,
		char const * const * const reads, char const * const * const refs,
		Alignment * const alignments) {

	int alignment_algorithm = opt & 0xF;

	fp_alignment_call alignment_call = 0;

	switch (alignment_algorithm) {
	case 0:
		alignment_call = &SSEKernel::calc_alignment_smith_waterman;
		//calc_alignment_smith_waterman(reads, refs, alignments);
		break;
	case 1:
		alignment_call = &SSEKernel::calc_alignment_needleman_wunsch;
		//calc_alignment_needleman_wunsch(reads, refs, alignments);
		break;
	default:
		// Unsupported mode
		break;
	}

	if (alignment_call != 0) {

		int num_batches = aln_number / SSE_SIZE;
		int mod = aln_number % SSE_SIZE;

		int cur_alignment = 0;

		Logger.log(0, KERNEL, "Started SSE aligning.");
		Logger.log(0, KERNEL,
				string("# Alignments:\t" + to_string(aln_number)).c_str());
		Logger.log(0, KERNEL,
				string("# Batches:\t" + to_string(num_batches)).c_str());
		Logger.log(0, KERNEL, string("# Overhang:\t" + to_string(mod)).c_str());

#pragma omp parallel for num_threads(Parameters.param_int("num_threads"))
		for (int i = num_batches; i > 0; --i) {
			(this->*alignment_call)(reads + cur_alignment, refs + cur_alignment,
					alignments + cur_alignment);
			cur_alignment += SSE_SIZE;
		}

		// Number of alignments not multiple of SSE_SIZE (8 alignments in SIMD registers)
		// -> Need to fillup remaining slots with \0 sequences

		if (mod != 0) {

			Alignment * alignment_overflow = new Alignment[SSE_SIZE];

			char const * * const read_overflow = new const char *[SSE_SIZE];
			char const * * const ref_overflow = new const char *[SSE_SIZE];

			char const * null_read = new char[readLength]();
			char const * null_ref = new char[refLength]();

			int overflow_fill = 0;

			while (overflow_fill < mod) {
				*(read_overflow + overflow_fill) = *(reads + cur_alignment
						+ overflow_fill);
				*(ref_overflow + overflow_fill) = *(refs + cur_alignment
						+ overflow_fill);
				++overflow_fill;
			}

			while (overflow_fill < SSE_SIZE) {
				*(read_overflow + overflow_fill) = null_read;
				*(ref_overflow + overflow_fill) = null_ref;
				++overflow_fill;
			}

			(this->*alignment_call)(read_overflow, ref_overflow,
					alignment_overflow);

			for (int i = 0; i < mod; ++i) {
				alignments[cur_alignment + i] = alignment_overflow[i];
			}

			delete[] alignment_overflow;
			alignment_overflow = 0;
			delete[] null_read;
			null_read = 0;
			delete[] null_ref;
			null_ref = 0;
			delete[] read_overflow;
			delete[] ref_overflow;
		}
	}
}

void SSEKernel::score_alignments(int const & opt, int const & aln_number,
		char const * const * const reads, char const * const * const refs,
		short * const scores) {

	/*###############################################################
	 * ##############################################################
	 * IF SEGFAULTS OCCUR TRY SAVING SCORES IN TEMPORY ALIGN16 ARRAY
	 * ##############################################################
	 * ##############################################################
	 */

	int alignment_algorithm = opt & 0xF;

	fp_scoring_call scoring_call = 0;

	switch (alignment_algorithm) {
	case 0:
		scoring_call = &SSEKernel::score_alignment_smith_waterman;
		break;
	case 1:
		scoring_call = &SSEKernel::score_alignment_needleman_wunsch;
		break;
	default:
		// Unsupported mode
		break;
	}

	if (scoring_call != 0) {

		int num_batches = aln_number / SSE_SIZE;
		int mod = aln_number % SSE_SIZE;

		int cur_alignment = 0;

		Logger.log(0, KERNEL, "Started SSE scoring.");
		Logger.log(0, KERNEL,
				string("# Alignments:\t" + to_string(aln_number)).c_str());
		Logger.log(0, KERNEL,
				string("# Batches:\t" + to_string(num_batches)).c_str());
		Logger.log(0, KERNEL, string("# Overhang:\t" + to_string(mod)).c_str());

		#pragma omp parallel for num_threads(Parameters.param_int("num_threads"))
		for (int i = num_batches; i > 0; --i) {
			(this->*scoring_call)(reads + cur_alignment, refs + cur_alignment,
					scores + cur_alignment);
			cur_alignment += SSE_SIZE;
		}

		// Number of alignments not multiple of SSE_SIZE (8 alignments in SIMD registers)
		// -> Need to fillup remaining slots with \0 sequences

		if (mod != 0) {

			short * score_overflow = new short[SSE_SIZE];

			char const * * const read_overflow = new const char *[SSE_SIZE];
			char const * * const ref_overflow = new const char *[SSE_SIZE];

			char const * null_read = new char[readLength]();
			char const * null_ref = new char[refLength]();

			int overflow_fill = 0;

			while (overflow_fill < mod) {
				*(read_overflow + overflow_fill) = *(reads + cur_alignment
						+ overflow_fill);
				*(ref_overflow + overflow_fill) = *(refs + cur_alignment
						+ overflow_fill);
				++overflow_fill;
			}

			while (overflow_fill < SSE_SIZE) {
				*(read_overflow + overflow_fill) = null_read;
				*(ref_overflow + overflow_fill) = null_ref;
				++overflow_fill;
			}
			(this->*scoring_call)(read_overflow, ref_overflow, score_overflow);

			for (int i = 0; i < mod; ++i) {
				scores[cur_alignment + i] = score_overflow[i];
			}

			delete[] score_overflow;
			score_overflow = 0;
			delete[] null_read;
			null_read = 0;
			delete[] null_ref;
			null_ref = 0;
			delete[] read_overflow;
			delete[] ref_overflow;
		}
	}
}

void SSEKernel::calc_alignment_matrix_smith_waterman(
		char const * const * const read, char const * const * const ref,
		short * const matrix, short * const best_coordinates) {

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

	malloc16(scoreMat,
			sizeof(short) * (refLength + 1) * SCORING_ROWS * SSE_SIZE, 16);
	memset(scoreMat, 0,
			sizeof(short) * (refLength + 1) * SCORING_ROWS * SSE_SIZE);

	// SSE conversion array for current read and ref bases
	align16 short read_bases[SSE_SIZE];
	align16 short ref_bases[SSE_SIZE];

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

		sse_read_bases = _mm_load_si128((__m128i   const *) read_bases);

		// UC read
		sse_read_bases = _mm_and_si128(sse_read_bases, x_UCMask);

#ifndef NDEBUG

		Logger.log(0, KERNEL,
				string(
						"Current read base:\t"
								+ __m128i_toString<char>(sse_read_bases)).c_str());

#endif

		__m128i valid_read_base = _mm_or_si128(
				_mm_cmpeq_epi16(sse_read_bases, x_C),
				_mm_or_si128(_mm_cmpeq_epi16(sse_read_bases, x_G),
						_mm_or_si128(_mm_cmpeq_epi16(sse_read_bases, x_T),
								_mm_cmpeq_epi16(sse_read_bases, x_A))));

		for (int ref_pos = 0; ref_pos < refLength; ++ref_pos) {

			// load ref base
			for (int base = 0; base < SSE_SIZE; ++base) {
				ref_bases[base] = *(ref[base] + ref_pos);
			}

			sse_ref_bases = _mm_load_si128((__m128i   const *) ref_bases);

			// load relevant matrix cells (up, diag, left)
			__m128i up = _mm_load_si128(
					(__m128i *) (scoreMat
							+ SSE_SIZE
									* (prev_row_score * (refLength + 1)
											+ ref_pos + 1)));
			__m128i diag = _mm_load_si128(
					(__m128i *) (scoreMat
							+ SSE_SIZE
									* (prev_row_score * (refLength + 1)
											+ ref_pos)));
			__m128i left = _mm_load_si128(
					(__m128i *) (scoreMat
							+ SSE_SIZE
									* (current_row_score * (refLength + 1)
											+ ref_pos)));

			// add gap penalties to up and left
			up = _mm_add_epi16(up, x_scoreGapRef);
			left = _mm_add_epi16(left, x_scoreGapRead);

			// UC ref
			sse_ref_bases = _mm_and_si128(sse_ref_bases, x_UCMask);

			__m128i valid_ref_base = _mm_or_si128(
					_mm_cmpeq_epi16(sse_ref_bases, x_C),
					_mm_or_si128(_mm_cmpeq_epi16(sse_ref_bases, x_G),
							_mm_or_si128(_mm_cmpeq_epi16(sse_ref_bases, x_T),
									_mm_cmpeq_epi16(sse_ref_bases, x_A))));

#ifndef NDEBUG

			Logger.log(0, KERNEL,
					string(
							"Current ref base:\t"
									+ __m128i_toString<char>(sse_ref_bases)).c_str());

#endif

			__m128i valid_comp = _mm_and_si128(valid_read_base, valid_ref_base);

			// match read and ref bases
			// matches will hold 1, mismatches will hold 0
			__m128i match = _mm_cmpeq_epi16(sse_read_bases, sse_ref_bases);

			// bitwise and between 0 (mismatches) and match score gives match score only for equal bases
			// add these to diagonal value

			diag = _mm_add_epi16(diag,
					_mm_and_si128(valid_comp,
							_mm_and_si128(match, x_scoreMatch)));
			//diag = _mm_add_epi16(diag, _mm_and_si128(match, x_scoreMatch));
			// bitwise not and between 1 (matches) and mismatch score gives mismatch score only for unequal bases
			// add these to diagonal value
			diag = _mm_add_epi16(diag,
					_mm_and_si128(valid_comp,
							_mm_andnot_si128(match, x_scoreMismatch)));

			// Cell value will be max of upper + gap penalty, left + gap penalty, diag + match/mismatch score or 0
			__m128i cell = _mm_max_epi16(diag,
					_mm_max_epi16(left, _mm_max_epi16(up, x_zeros)));

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
			match = _mm_and_si128(valid_comp, _mm_cmpeq_epi16(cell, diag));
			pointer = _mm_max_epi16(pointer, _mm_and_si128(match, x_p_diag));

			// Store score and pointer
			_mm_store_si128(
					(__m128i *) (scoreMat
							+ SSE_SIZE
									* (current_row_score * (refLength + 1)
											+ ref_pos + 1)), cell);
			_mm_store_si128(
					(__m128i *) (matrix
							+ SSE_SIZE
									* (current_row_aln * (refLength + 1)
											+ ref_pos + 1)), pointer);

			// Store read and ref positions if new max score
			match = _mm_cmpgt_epi16(cell, max_score);

#ifndef NDEBUG

			Logger.log(0, KERNEL,
					string("Scores:\t" + __m128i_toString<short>(cell)).c_str());
			Logger.log(0, KERNEL,
					string("Max:\t" + __m128i_toString<short>(max_score)).c_str());
			Logger.log(0, KERNEL,
					string("Match:\t" + __m128i_toString<short>(match)).c_str());

#endif

			best_read_pos = _mm_max_epi16(
					_mm_andnot_si128(match, best_read_pos),
					_mm_and_si128(match, sse_read_pos));
			best_ref_pos = _mm_max_epi16(_mm_andnot_si128(match, best_ref_pos),
					_mm_and_si128(match, sse_ref_pos));

#ifndef NDEBUG

			Logger.log(0, KERNEL,
					string(
							"Best read pos:\t"
									+ __m128i_toString<short>(best_read_pos)).c_str());
			Logger.log(0, KERNEL,
					string(
							"Best ref pos:\t"
									+ __m128i_toString<short>(best_ref_pos)).c_str());

#endif

			max_score = _mm_max_epi16(cell, max_score);

			sse_ref_pos = _mm_add_epi16(sse_ref_pos, see_increment);

		}
		prev_row_score = current_row_score;
		(++current_row_score) &= 1;

		++current_row_aln;

		sse_read_pos = _mm_add_epi16(sse_read_pos, see_increment);

	}

#ifndef NDEBUG

	Logger.log(0, KERNEL,
			string("Max scores:\t" + __m128i_toString<short>(max_score)).c_str());
	Logger.log(0, KERNEL,
			string("Best ref pos:\t" + __m128i_toString<short>(best_ref_pos)).c_str());

#endif

	free(scoreMat);

	_mm_store_si128((__m128i *) best_coordinates, best_read_pos);
	_mm_store_si128((__m128i *) best_coordinates + 1, best_ref_pos);
}

void SSEKernel::calculate_alignment_matrix_needleman_wunsch(
		char const * const * const read, char const * const * const ref,
		short * const matrix, short * const best_coordinates) {

	// Tracking positions where invalid characters start
	__m128i max_read_pos = short_to_sse(readLength - 1);
	__m128i max_ref_pos = short_to_sse(refLength - 1);

	// Tracking rowwise maxima
	__m128i rowMax = short_to_sse(SHRT_MIN);
	__m128i rowMaxIndex = x_zeros;

	// Tracking global maxima
	//	__m128i globalRowMax = short_to_sse(SHRT_MIN);
	__m128i globalRowMaxIndex = short_to_sse(-1);

	// Scoring matrix indices
	short prev_row_score = 0;
	short current_row_score = 1;

	// Alignment matrix index
	short current_row_aln = 1;

	short * scoreMat = 0;

	malloc16(scoreMat,
			sizeof(short) * (refLength + 1) * SCORING_ROWS * SSE_SIZE, 16);
	memset(scoreMat, 0,
			sizeof(short) * (refLength + 1) * SCORING_ROWS * SSE_SIZE);

	// SSE conversion array for current read and ref bases
	align16 short read_bases[SSE_SIZE];
	align16 short ref_bases[SSE_SIZE];

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

		sse_read_bases = _mm_load_si128((__m128i   const *) read_bases);

		// UC read
		sse_read_bases = _mm_and_si128(sse_read_bases, x_UCMask);

		__m128i valid_read_base = _mm_or_si128(
				_mm_cmpeq_epi16(sse_read_bases, x_C),
				_mm_or_si128(_mm_cmpeq_epi16(sse_read_bases, x_G),
						_mm_or_si128(_mm_cmpeq_epi16(sse_read_bases, x_T),
								_mm_cmpeq_epi16(sse_read_bases, x_A))));

		// First alignment column is always UP
		_mm_store_si128(
				(__m128i *) (matrix
						+ SSE_SIZE * (current_row_aln * (refLength + 1))),
				x_p_up);
		// First score column is always continued gap-ref score
		_mm_store_si128(
				(__m128i *) (scoreMat
						+ SSE_SIZE * (current_row_score * (refLength + 1))),
				short_to_sse((read_pos + 1) * scoreGapRef));

		// if max_readpos == readLength - 1 AND read char is invalid -> max_read_pos = read_pos - 1
		max_read_pos = _mm_blendv_si128(max_read_pos,
				_mm_sub_epi16(sse_read_pos, see_increment),
				_mm_andnot_si128(valid_read_base,
						_mm_cmpeq_epi16(max_read_pos,
								short_to_sse(readLength - 1))));

#ifndef NDEBUG

		Logger.log(0, KERNEL,
				string(
						"Current read base:\t"
								+ __m128i_toString<char>(sse_read_bases)).c_str());

		Logger.log(0, KERNEL,
				string(
						"Max read pos:\t"
								+ __m128i_toString<short>(max_read_pos)).c_str());
		Logger.log(0, KERNEL,
				string("Read pos:\t" + __m128i_toString<short>(sse_read_pos)).c_str());

#endif

		// Save previous row max if read ends prematurely
		// if max_readpos + 1 == read_pos -> globalRowMax = rowMax; globalRowMaxIndex = rowMaxIndex;
		globalRowMaxIndex = _mm_blendv_si128(globalRowMaxIndex, rowMaxIndex,
				_mm_cmpeq_epi16(_mm_add_epi16(max_read_pos, see_increment),
						sse_read_pos));

#ifndef NDEBUG

		Logger.log(0, KERNEL,
				string(
						"global_row_max_index:\t"
								+ __m128i_toString<short>(globalRowMaxIndex)).c_str());

#endif

		rowMax = _mm_load_si128(
				(__m128i *) (scoreMat
						+ SSE_SIZE * (current_row_score * (refLength + 1))));
		rowMaxIndex = x_zeros;

		for (int ref_pos = 0; ref_pos < refLength; ++ref_pos) {

			// load ref base
			for (int base = 0; base < SSE_SIZE; ++base) {
				ref_bases[base] = *(ref[base] + ref_pos);
			}

			sse_ref_bases = _mm_load_si128((__m128i   const *) ref_bases);

			// load relevant matrix cells (up, diag, left)
			__m128i up = _mm_load_si128(
					(__m128i *) (scoreMat
							+ SSE_SIZE
									* (prev_row_score * (refLength + 1)
											+ ref_pos + 1)));
			__m128i diag = _mm_load_si128(
					(__m128i *) (scoreMat
							+ SSE_SIZE
									* (prev_row_score * (refLength + 1)
											+ ref_pos)));
			__m128i left = _mm_load_si128(
					(__m128i *) (scoreMat
							+ SSE_SIZE
									* (current_row_score * (refLength + 1)
											+ ref_pos)));

			// add gap penalties to up and left
			up = _mm_add_epi16(up, x_scoreGapRef);
			left = _mm_add_epi16(left, x_scoreGapRead);

			// UC ref
			sse_ref_bases = _mm_and_si128(sse_ref_bases, x_UCMask);

			__m128i valid_ref_base = _mm_or_si128(
					_mm_cmpeq_epi16(sse_ref_bases, x_C),
					_mm_or_si128(_mm_cmpeq_epi16(sse_ref_bases, x_G),
							_mm_or_si128(_mm_cmpeq_epi16(sse_ref_bases, x_T),
									_mm_cmpeq_epi16(sse_ref_bases, x_A))));

#ifndef NDEBUG

			Logger.log(0, KERNEL,
					string(
							"Current ref base:\t"
									+ __m128i_toString<char>(sse_ref_bases)).c_str());

#endif

			__m128i valid_comp = _mm_and_si128(valid_read_base, valid_ref_base);

			// match read and ref bases
			// matches will hold 1, mismatches will hold 0
			__m128i match = _mm_cmpeq_epi16(sse_read_bases, sse_ref_bases);

			// bitwise and between 0 (mismatches) and match score gives match score only for equal bases
			// add these to diagonal value

			diag = _mm_add_epi16(diag,
					_mm_and_si128(valid_comp,
							_mm_and_si128(match, x_scoreMatch)));
			//diag = _mm_add_epi16(diag, _mm_and_si128(match, x_scoreMatch));
			// bitwise not and between 1 (matches) and mismatch score gives mismatch score only for unequal bases
			// add these to diagonal value
			diag = _mm_add_epi16(diag,
					_mm_and_si128(valid_comp,
							_mm_andnot_si128(match, x_scoreMismatch)));

			// Cell value will be max of upper + gap penalty, left + gap penalty, diag + match/mismatch score or 0
			__m128i cell = _mm_max_epi16(diag, _mm_max_epi16(left, up));

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
			match = _mm_and_si128(valid_comp, _mm_cmpeq_epi16(cell, diag));
			pointer = _mm_max_epi16(pointer, _mm_and_si128(match, x_p_diag));

			// Store score and pointer
			_mm_store_si128(
					(__m128i *) (scoreMat
							+ SSE_SIZE
									* (current_row_score * (refLength + 1)
											+ ref_pos + 1)), cell);
			_mm_store_si128(
					(__m128i *) (matrix
							+ SSE_SIZE
									* (current_row_aln * (refLength + 1)
											+ ref_pos + 1)), pointer);

			max_ref_pos = _mm_blendv_si128(max_ref_pos,
					_mm_sub_epi16(sse_ref_pos, see_increment),
					_mm_andnot_si128(valid_ref_base,
							_mm_cmpeq_epi16(max_ref_pos,
									short_to_sse(refLength - 1))));

#ifndef NDEBUG

			Logger.log(0, KERNEL,
					string("Scores:\t" + __m128i_toString<short>(cell)).c_str());
			Logger.log(0, KERNEL,
					string(
							"Max_ref_pos:\t"
									+ __m128i_toString<short>(max_ref_pos)).c_str());

#endif

			// if cur > rowMax => rowMax = cur; rowMaxIndex = ref_pos;
			__m128i sel = _mm_cmpgt_epi16(cell, rowMax);
			rowMax = _mm_blendv_si128(rowMax, cell, sel);
			rowMaxIndex = _mm_blendv_si128(rowMaxIndex, sse_ref_pos, sel);

			sse_ref_pos = _mm_add_epi16(sse_ref_pos, see_increment);

		}
		prev_row_score = current_row_score;
		(++current_row_score) &= 1;

		++current_row_aln;

		sse_read_pos = _mm_add_epi16(sse_read_pos, see_increment);

	}

	free(scoreMat);

	globalRowMaxIndex = _mm_blendv_si128(globalRowMaxIndex, rowMaxIndex,
			_mm_cmplt_epi16(globalRowMaxIndex, x_zeros));

#ifndef NDEBUG

	Logger.log(0, KERNEL,
			string("Max_ref_pos:\t" + __m128i_toString<short>(max_ref_pos)).c_str());
	Logger.log(0, KERNEL,
			string(
					"global_row_max_index:\t"
							+ __m128i_toString<short>(globalRowMaxIndex)).c_str());

#endif

	__m128i best_ref_pos = _mm_min_epi16(max_ref_pos, globalRowMaxIndex);

	_mm_store_si128((__m128i *) best_coordinates, max_read_pos);
	_mm_store_si128((__m128i *) best_coordinates + 1, best_ref_pos);
}

void SSEKernel::calc_alignment_smith_waterman(char const * const * const read,
		char const * const * const ref, Alignment * const alignment) {

	short * matrix = 0;
	malloc16(matrix,
			sizeof(short) * (refLength + 1) * (readLength + 1) * SSE_SIZE, 16);
	memset(matrix, 0,
			(refLength + 1) * (readLength + 1) * SSE_SIZE * sizeof(short));

	short * best_coordinates = 0;
	malloc16(best_coordinates, sizeof(short) * 2 * SSE_SIZE, 16);
	memset(best_coordinates, 0, sizeof(short) * 2 * SSE_SIZE);

	calc_alignment_matrix_smith_waterman(read, ref, matrix, best_coordinates);


#ifndef NDEBUG

	Logger.log(0, KERNEL, "Score matrix:");

	for (int SSE_register = 0; SSE_register < SSE_SIZE; ++SSE_register) {

		Logger.log(0, KERNEL,
				string("Alignment matrix " + to_string(SSE_register) + ":").c_str());
		for (int i = 0; i < readLength + 1; ++i) {
			string matrix_line;
			for (int j = 0; j < refLength + 1; ++j) {
				matrix_line += (to_string(
						*(matrix + SSE_SIZE * (i * (refLength + 1) + j)
								+ SSE_register)) + "\t");
			}
			Logger.log(0, KERNEL, matrix_line.c_str());
		}
		Logger.log(0, KERNEL, "");
	}

#endif

	char * alignments = new char[alnLength * 2 * SSE_SIZE];

	// Retreive best read and ref pos
	__m128i best_read_positions = _mm_load_si128((__m128i *) best_coordinates);
	__m128i best_ref_positions = _mm_load_si128(
			(__m128i *) best_coordinates + 1);

#ifndef NDEBUG

	Logger.log(0, KERNEL,
			string(
					"Best read:\t"
							+ __m128i_toString<short>(best_read_positions)).c_str());
	Logger.log(0, KERNEL,
			string("Best ref:\t" + __m128i_toString<short>(best_ref_positions)).c_str());

#endif

	for (int SSE_register = 0; SSE_register < SSE_SIZE; ++SSE_register) {

		short read_pos = *(best_coordinates + SSE_register);
		short ref_pos = *(best_coordinates + SSE_SIZE + SSE_register);

		int aln_pos = alnLength - 2;

		alignments[(SSE_register * alnLength * 2) + alnLength - 1] = '\0';
		alignments[(SSE_register * alnLength * 2) + (2 * alnLength) - 1] = '\0';

		short backtrack = *(matrix
				+ SSE_SIZE * ((read_pos + 1) * (refLength + 1) + ref_pos + 1)
				+ SSE_register);

		while (backtrack != START) {

#ifndef NDEBUG

			Logger.log(0, KERNEL,
					string("Current backtrack:\t" + to_string(backtrack)).c_str());

#endif

			if (backtrack == UP) {
				alignments[(SSE_register * alnLength * 2) + alnLength + aln_pos] =
						'-';
				alignments[(SSE_register * alnLength * 2) + aln_pos] =
						read[SSE_register][read_pos--];
			}

			if (backtrack == LEFT) {
				alignments[(SSE_register * alnLength * 2) + aln_pos] = '-';
				alignments[(SSE_register * alnLength * 2) + alnLength + aln_pos] =
						ref[SSE_register][ref_pos--];
			}

			if (backtrack == DIAG) {
				alignments[(SSE_register * alnLength * 2) + aln_pos] =
						read[SSE_register][read_pos--];
				alignments[(SSE_register * alnLength * 2) + alnLength + aln_pos] =
						ref[SSE_register][ref_pos--];
			}

			backtrack = *(matrix
					+ SSE_SIZE
							* ((read_pos + 1) * (refLength + 1) + ref_pos + 1)
					+ SSE_register);
			--aln_pos;

#ifndef NDEBUG

			Logger.log(0, KERNEL,
					string("Current read pos:\t" + to_string(read_pos)).c_str());
			Logger.log(0, KERNEL,
					string("Current ref pos:\t" + to_string(ref_pos)).c_str());

#endif

		}

		alignment[SSE_register].read = new char[alnLength];
		alignment[SSE_register].ref = new char[alnLength];

		memcpy(alignment[SSE_register].read,
				alignments + (SSE_register * alnLength * 2),
				alnLength * sizeof(char));
		memcpy(alignment[SSE_register].ref,
				alignments + (SSE_register * alnLength * 2) + alnLength,
				alnLength * sizeof(char));

		alignment[SSE_register].readStart = aln_pos + 1;
		alignment[SSE_register].refStart = aln_pos + 1;

		alignment[SSE_register].readEnd = alnLength - 1;
		alignment[SSE_register].refEnd = alnLength - 1;
	}

	delete[] alignments;
	alignments = 0;
	free(best_coordinates);
	free(matrix);
}

void SSEKernel::calc_alignment_needleman_wunsch(char const * const * const read,
		char const * const * const ref, Alignment * const alignment) {

	short * matrix = 0;
	malloc16(matrix,
			sizeof(short) * (refLength + 1) * (readLength + 1) * SSE_SIZE, 16);
	memset(matrix, 0,
			(refLength + 1) * (readLength + 1) * SSE_SIZE * sizeof(short));

	short * best_coordinates = 0;
	malloc16(best_coordinates, sizeof(short) * 2 * SSE_SIZE, 16);
	memset(best_coordinates, 0, sizeof(short) * 2 * SSE_SIZE);

	calculate_alignment_matrix_needleman_wunsch(read, ref, matrix,
			best_coordinates);

#ifndef NDEBUG

	Logger.log(0, KERNEL, "Score matrix:");

	for (int SSE_register = 0; SSE_register < SSE_SIZE; ++SSE_register) {

		Logger.log(0, KERNEL,
				string("Alignment matrix " + to_string(SSE_register) + ":").c_str());
		for (int i = 0; i < readLength + 1; ++i) {
			string matrix_line;
			for (int j = 0; j < refLength + 1; ++j) {
				matrix_line += (to_string(
						*(matrix + SSE_SIZE * (i * (refLength + 1) + j)
								+ SSE_register)) + "\t");
			}
			Logger.log(0, KERNEL, matrix_line.c_str());
		}
		Logger.log(0, KERNEL, "");
	}

#endif

	char * alignments = new char[alnLength * 2 * SSE_SIZE];

	// Retreive best read and ref pos
	__m128i best_read_positions = _mm_load_si128((__m128i *) best_coordinates);
	__m128i best_ref_positions = _mm_load_si128(
			(__m128i *) best_coordinates + 1);

#ifndef NDEBUG

	Logger.log(0, KERNEL,
			string(
					"Best read:\t"
							+ __m128i_toString<short>(best_read_positions)).c_str());
	Logger.log(0, KERNEL,
			string("Best ref:\t" + __m128i_toString<short>(best_ref_positions)).c_str());

#endif

	for (int SSE_register = 0; SSE_register < SSE_SIZE; ++SSE_register) {

		short read_pos = *(best_coordinates + SSE_register);
		short ref_pos = *(best_coordinates + SSE_SIZE + SSE_register);

		int aln_pos = alnLength - 2;

		alignments[(SSE_register * alnLength * 2) + alnLength - 1] = '\0';
		alignments[(SSE_register * alnLength * 2) + (2 * alnLength) - 1] = '\0';

		short backtrack = *(matrix
				+ SSE_SIZE * ((read_pos + 1) * (refLength + 1) + ref_pos + 1)
				+ SSE_register);

		while (backtrack != START) {

#ifndef NDEBUG

			Logger.log(0, KERNEL,
					string("Current backtrack:\t" + to_string(backtrack)).c_str());

#endif

			if (backtrack == UP) {
				alignments[(SSE_register * alnLength * 2) + alnLength + aln_pos] =
						'-';
				alignments[(SSE_register * alnLength * 2) + aln_pos] =
						read[SSE_register][read_pos--];
			}

			if (backtrack == LEFT) {
				alignments[(SSE_register * alnLength * 2) + aln_pos] = '-';
				alignments[(SSE_register * alnLength * 2) + alnLength + aln_pos] =
						ref[SSE_register][ref_pos--];
			}

			if (backtrack == DIAG) {
				alignments[(SSE_register * alnLength * 2) + aln_pos] =
						read[SSE_register][read_pos--];
				alignments[(SSE_register * alnLength * 2) + alnLength + aln_pos] =
						ref[SSE_register][ref_pos--];
			}

			backtrack = *(matrix
					+ SSE_SIZE
							* ((read_pos + 1) * (refLength + 1) + ref_pos + 1)
					+ SSE_register);
			--aln_pos;

#ifndef NDEBUG

			Logger.log(0, KERNEL,
					string("Current read pos:\t" + to_string(read_pos)).c_str());
			Logger.log(0, KERNEL,
					string("Current ref pos:\t" + to_string(ref_pos)).c_str());

#endif
		}

		alignment[SSE_register].read = new char[alnLength];
		alignment[SSE_register].ref = new char[alnLength];

		memcpy(alignment[SSE_register].read,
				alignments + (SSE_register * alnLength * 2),
				alnLength * sizeof(char));
		memcpy(alignment[SSE_register].ref,
				alignments + (SSE_register * alnLength * 2) + alnLength,
				alnLength * sizeof(char));

		alignment[SSE_register].readStart = aln_pos + 1;
		alignment[SSE_register].refStart = aln_pos + 1;

		alignment[SSE_register].readEnd = alnLength - 1;
		alignment[SSE_register].refEnd = alnLength - 1;
	}

	delete[] alignments;
	alignments = 0;
	free(best_coordinates);
	free(matrix);

}

void SSEKernel::score_alignment_smith_waterman(char const * const * const read,
		char const * const * const ref, short * const scores) {

	__m128i max_score = x_zeros;

	// Initialize SSE matrix
	short * matrix = 0;

	malloc16(matrix, sizeof(short) * (refLength + 1) * SCORING_ROWS * SSE_SIZE,
			16);
	memset(matrix, 0,
			sizeof(short) * (refLength + 1) * SCORING_ROWS * SSE_SIZE);

	// offset for first and second row
	int prev_row = 0;
	int cur_row = 1;

	// SSE conversion array for current read and ref bases
	align16 short read_bases[SSE_SIZE];
	align16 short ref_bases[SSE_SIZE];

	// holds current read and ref base for SSE instruction
	__m128i sse_read_bases;
	__m128i sse_ref_bases;

	for (int read_pos = 0; read_pos < readLength; ++read_pos) {

		// load read base
		for (int base = 0; base < SSE_SIZE; ++base) {
			read_bases[base] = *(read[base] + read_pos);
		}

		sse_read_bases = _mm_load_si128((__m128i   const *) read_bases);

		// UC read
		sse_read_bases = _mm_and_si128(sse_read_bases, x_UCMask);

#ifndef NDEBUG

			Logger.log(0, KERNEL,string("Read base:\t" + __m128i_toString<char>(sse_read_bases)).c_str());

#endif

		__m128i valid_read_base = _mm_or_si128(
				_mm_cmpeq_epi16(sse_read_bases, x_C),
				_mm_or_si128(_mm_cmpeq_epi16(sse_read_bases, x_G),
						_mm_or_si128(_mm_cmpeq_epi16(sse_read_bases, x_T),
								_mm_cmpeq_epi16(sse_read_bases, x_A))));

		for (int ref_pos = 0; ref_pos < refLength; ++ref_pos) {

			// load ref base
			for (int base = 0; base < SSE_SIZE; ++base) {
				ref_bases[base] = *(ref[base] + ref_pos);
			}

			sse_ref_bases = _mm_load_si128((__m128i   const *) ref_bases);

			// load relevant matrix cells (up, diag, left)
			__m128i up =
					_mm_load_si128(
							(__m128i *) (matrix
									+ SSE_SIZE
											* (prev_row * (refLength + 1)
													+ ref_pos + 1)));
			__m128i diag =
					_mm_load_si128(
							(__m128i *) (matrix
									+ SSE_SIZE
											* (prev_row * (refLength + 1)
													+ ref_pos)));
			__m128i left =
					_mm_load_si128(
							(__m128i *) (matrix
									+ SSE_SIZE
											* (cur_row * (refLength + 1)
													+ ref_pos)));

			// add gap penalties to up and left
			up = _mm_add_epi16(up, x_scoreGapRef);
			left = _mm_add_epi16(left, x_scoreGapRead);

			// UC ref
			sse_ref_bases = _mm_and_si128(sse_ref_bases, x_UCMask);

			__m128i valid_ref_base = _mm_or_si128(
					_mm_cmpeq_epi16(sse_ref_bases, x_C),
					_mm_or_si128(_mm_cmpeq_epi16(sse_ref_bases, x_G),
							_mm_or_si128(_mm_cmpeq_epi16(sse_ref_bases, x_T),
									_mm_cmpeq_epi16(sse_ref_bases, x_A))));

			__m128i valid_comp = _mm_and_si128(valid_read_base, valid_ref_base);

#ifndef NDEBUG

			Logger.log(0, KERNEL,string("Ref base:\t" + __m128i_toString<char>(sse_ref_bases)).c_str());

#endif

			// match read and ref bases
			// matches will hold 1, mismatches will hold 0
			__m128i match = _mm_cmpeq_epi16(sse_read_bases, sse_ref_bases);

			// bitwise and between 0 (mismatches) and match score gives match score only for equal bases
			// add these to diagonal value
			diag = _mm_add_epi16(diag,
					_mm_and_si128(valid_comp,
							_mm_and_si128(match, x_scoreMatch)));
			//diag = _mm_add_epi16(diag, _mm_and_si128(match, x_scoreMatch));
			// bitwise not and between 1 (matches) and mismatch score gives mismatch score only for unequal bases
			// add these to diagonal value
			//diag = _mm_add_epi16(diag, _mm_andnot_si128(match, x_scoreMismatch));
			diag = _mm_add_epi16(diag,
					_mm_and_si128(valid_comp,
							_mm_andnot_si128(match, x_scoreMismatch)));

			// Cell value will be max of upper + gap penalty, left + gap penalty, diag + match/mismatch score or 0
			__m128i cell = _mm_max_epi16(diag,
					_mm_max_epi16(left, _mm_max_epi16(up, x_zeros)));

			_mm_store_si128(
					(__m128i *) (matrix
							+ SSE_SIZE
									* (cur_row * (refLength + 1) + ref_pos + 1)),
					cell);

			max_score = _mm_max_epi16(cell, max_score);

#ifndef NDEBUG

			Logger.log(0, KERNEL,string("Cell:\t" + __m128i_toString<short>(cell)).c_str());
			Logger.log(0, KERNEL,string("Max:\t" + __m128i_toString<short>(max_score)).c_str());

#endif

		}
		prev_row = cur_row;
		(++cur_row) &= 1;
	}

	free(matrix);

	_mm_store_si128((__m128i *) scores, max_score);
}

void SSEKernel::score_alignment_needleman_wunsch(
		char const * const * const read, char const * const * const ref,
		short * const scores) {

	__m128i max_score = x_zeros;

	// Initialize SSE matrix
	short * matrix = 0;

	malloc16(matrix, sizeof(short) * (refLength + 1) * SCORING_ROWS * SSE_SIZE,
			16);
	memset(matrix, 0,
			sizeof(short) * (refLength + 1) * SCORING_ROWS * SSE_SIZE);

	// offset for first and second row
	int prev_row = 0;
	int cur_row = 1;

	// SSE conversion array for current read and ref bases
	align16 short read_bases[SSE_SIZE];
	align16 short ref_bases[SSE_SIZE];

	// holds current read and ref base for SSE instruction
	__m128i sse_read_bases;
	__m128i sse_ref_bases;

	for (int read_pos = 0; read_pos < readLength; ++read_pos) {

		// load read base
		for (int base = 0; base < SSE_SIZE; ++base) {
			read_bases[base] = *(read[base] + read_pos);
		}

		sse_read_bases = _mm_load_si128((__m128i   const *) read_bases);

		// UC read
		sse_read_bases = _mm_and_si128(sse_read_bases, x_UCMask);

#ifndef NDEBUG

			Logger.log(0, KERNEL,string("Read base:\t" + __m128i_toString<char>(sse_read_bases)).c_str());

#endif

		__m128i valid_read_base = _mm_or_si128(
				_mm_cmpeq_epi16(sse_read_bases, x_C),
				_mm_or_si128(_mm_cmpeq_epi16(sse_read_bases, x_G),
						_mm_or_si128(_mm_cmpeq_epi16(sse_read_bases, x_T),
								_mm_cmpeq_epi16(sse_read_bases, x_A))));

		for (int ref_pos = 0; ref_pos < refLength; ++ref_pos) {

			// load ref base
			for (int base = 0; base < SSE_SIZE; ++base) {
				ref_bases[base] = *(ref[base] + ref_pos);
			}

			sse_ref_bases = _mm_load_si128((__m128i   const *) ref_bases);

			// load relevant matrix cells (up, diag, left)
			__m128i up =
					_mm_load_si128(
							(__m128i *) (matrix
									+ SSE_SIZE
											* (prev_row * (refLength + 1)
													+ ref_pos + 1)));
			__m128i diag =
					_mm_load_si128(
							(__m128i *) (matrix
									+ SSE_SIZE
											* (prev_row * (refLength + 1)
													+ ref_pos)));
			__m128i left =
					_mm_load_si128(
							(__m128i *) (matrix
									+ SSE_SIZE
											* (cur_row * (refLength + 1)
													+ ref_pos)));

			// add gap penalties to up and left
			up = _mm_add_epi16(up, x_scoreGapRef);
			left = _mm_add_epi16(left, x_scoreGapRead);

			// UC ref
			sse_ref_bases = _mm_and_si128(sse_ref_bases, x_UCMask);

			__m128i valid_ref_base = _mm_or_si128(
					_mm_cmpeq_epi16(sse_ref_bases, x_C),
					_mm_or_si128(_mm_cmpeq_epi16(sse_ref_bases, x_G),
							_mm_or_si128(_mm_cmpeq_epi16(sse_ref_bases, x_T),
									_mm_cmpeq_epi16(sse_ref_bases, x_A))));

			__m128i valid_comp = _mm_and_si128(valid_read_base, valid_ref_base);

#ifndef NDEBUG

			Logger.log(0, KERNEL,string("Ref base:\t" + __m128i_toString<char>(sse_ref_bases)).c_str());

#endif

			// match read and ref bases
			// matches will hold 1, mismatches will hold 0
			__m128i match = _mm_cmpeq_epi16(sse_read_bases, sse_ref_bases);

			// bitwise and between 0 (mismatches) and match score gives match score only for equal bases
			// add these to diagonal value
			diag = _mm_add_epi16(diag,
					_mm_and_si128(valid_comp,
							_mm_and_si128(match, x_scoreMatch)));
			//diag = _mm_add_epi16(diag, _mm_and_si128(match, x_scoreMatch));
			// bitwise not and between 1 (matches) and mismatch score gives mismatch score only for unequal bases
			// add these to diagonal value
			//diag = _mm_add_epi16(diag, _mm_andnot_si128(match, x_scoreMismatch));
			diag = _mm_add_epi16(diag,
					_mm_and_si128(valid_comp,
							_mm_andnot_si128(match, x_scoreMismatch)));

			// Cell value will be max of upper + gap penalty, left + gap penalty, diag + match/mismatch score or 0
			__m128i cell = _mm_max_epi16(diag, _mm_max_epi16(left, up));

			_mm_store_si128(
					(__m128i *) (matrix
							+ SSE_SIZE
									* (cur_row * (refLength + 1) + ref_pos + 1)),
					cell);

#ifndef NDEBUG

			Logger.log(0, KERNEL,string("Cell:\t" + __m128i_toString<short>(cell)).c_str());

#endif
		}

		max_score =
				_mm_max_epi16(max_score,
						_mm_load_si128(
								(__m128i *) (matrix
										+ SSE_SIZE
												* (cur_row * (refLength + 1)
														+ refLength))));
#ifndef NDEBUG

		Logger.log(0, KERNEL,string("Cell:\t" + __m128i_toString<short>(max_score)).c_str());

#endif

		prev_row = cur_row;
		(++cur_row) &= 1;
	}

	for (int ref_pos = 0; ref_pos < refLength + 1; ++ref_pos) {
		max_score =
				_mm_max_epi16(max_score,
						_mm_load_si128(
								(__m128i *) (matrix
										+ SSE_SIZE
												* (prev_row * (refLength + 1)
														+ ref_pos))));
	}

	free(matrix);

	_mm_store_si128((__m128i *) scores, max_score);
}

#undef SCORING_ROWS
#undef KERNEL
#undef SSE_SIZE
