/*
 * AVXKernel.h
 *
 *  Created on: May 23, 2016
 *      Author: tobias.neumann
 */

#ifndef AVXKERNEL_H
#define AVXKERNEL_H

#include "AlignmentKernel.h"
#include "AlignmentParameters.h"

#if _WIN32
#define align32 __declspec(align(32))
#define malloc32(ptr,size)  ptr = (short*) _aligned_malloc(size, 32)
#else
#include <stdlib.h>
#define malloc32(ptr,size)  posix_memalign(((void * *)&ptr), 32, size)
#define align32 __attribute__((aligned(32)))
#endif

#define AVX_SIZE 16

#define UP 1
#define LEFT 2
#define DIAG 3
#define START 0

#include <immintrin.h>
#include <climits>

#include <iostream>

class AVXKernel : public AlignmentKernel{

public:
	AVXKernel() {

		bool exception = false;

		Parameters.has_key("score_match") ? scoreMatch = Parameters.param_int("score_match") : exception = true;
		Parameters.has_key("score_mismatch") ? scoreMismatch = Parameters.param_int("score_mismatch") : exception = true;
		Parameters.has_key("score_gap_read") ? scoreGapRead = Parameters.param_int("score_gap_read") : exception = true;
		Parameters.has_key("score_gap_ref") ? scoreGapRef = Parameters.param_int("score_gap_ref") : exception = true;
		Parameters.has_key("read_length") ? readLength = Parameters.param_int("read_length") : exception = true;
		Parameters.has_key("ref_length") ? refLength = Parameters.param_int("ref_length") : exception = true;

		alnLength = refLength + readLength;

		std::cout << "Match: " << scoreMatch
						<< "\nMismatch: " << scoreMismatch
						<< "\nGap_read: " << scoreGapRead
						<< "\nGap_ref: " << scoreGapRef
						<< "\nRead_length: " << readLength
						<< "\nRef_length: " << refLength
						<< "\nAln_length: " << alnLength << std::endl;

		if (exception) {
			throw "Cannot instantiate Kernel. Lacking parameters";
		}

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

	virtual ~AVXKernel() {}

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

	void calc_alignment_needleman_wunsch(char const * const * const read,
			char const * const * const ref, Alignment * const alignment);

	void calc_alignment_matrix_smith_waterman(char const * const * const read,
			char const * const * const ref, short * const matrix, short * const best_coordinates);

	void calculate_alignment_matrix_needleman_wunsch(char const * const * const read,
			char const * const * const ref, short * const matrix, short * const best_coordinates);

	typedef void (AVXKernel::* fp_alignment_call)(char const * const * const,  char const * const * const, Alignment * const);
	typedef void (AVXKernel::* fp_scoring_call)(char const * const * const,  char const * const * const, short * const);

	// Short = 2 byte
	// __m256i fits 256 bits = 16 shorts
	inline __m256i short_to_avx(short x) {

		short buf[AVX_SIZE] align32;
		for (int i = 0; i < AVX_SIZE; ++i) {
			buf[i] = x;
//			std::cout << "Saving " << x << "\n";
//					std::cin >> a;
		}

		__m256i ret = _mm256_load_si256((__m256i *) buf);
		return ret;
	}

	inline __m256i _mm_blendv_si256 (__m256i x, __m256i y, __m256i mask) {
		// Replace bit in x with bit in y when matching bit in mask is set:
		return _mm256_or_si256(_mm256_andnot_si256(mask, x), _mm256_and_si256(mask, y));
	}

	int readLength;
	int refLength;
	int alnLength;

	align32 short scoreMatch;
	align32 short scoreMismatch;
	align32 short scoreGapRead;
	align32 short scoreGapRef;

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

#endif /* AVXKERNEL_H */
