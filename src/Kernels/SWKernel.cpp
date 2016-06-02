/*
 * SWAligner.cpp
 *
 *  Created on: May 23, 2016
 *      Author: tobias.neumann
 */

#include "SWKernel.h"

using std::cout;
using std::endl;

inline short max(short a, short b) {
  return a > b ? a : b;
}

// Class method implementation

SWKernel::~SWKernel() {
	// TODO Auto-generated destructor stub
}

// Main functions

void SWKernel::score_alignment(char const * const * const read,
		char const * const * const ref, short * const scores) {

	short max_score = 0;

	short * matrix = new short[(refLength + 1) * 2]();

	int prev_row = 0;
	int current_row = 1;

	for (int read_pos = 0; read_pos < readLength; ++read_pos) {

			for (int ref_pos = 0; ref_pos < refLength; ++ref_pos) {

				short up = matrix[prev_row * (refLength + 1) + ref_pos + 1];
				short diag = matrix[prev_row * (refLength + 1) + ref_pos];
				short left = matrix[current_row * (refLength + 1) + ref_pos];

				read[0][read_pos] == ref[0][ref_pos] ? diag += scoreMatch : diag += scoreMismatch;

				short cur = max(up + scoreGapRef, max(left + scoreGapRead, max(diag, 0)));

				matrix[current_row * (refLength + 1) + ref_pos + 1] = cur;

				max_score = max(max_score, cur);

			}
			prev_row = current_row;
			(++current_row) &= 1;
		}
	delete [] matrix; matrix = 0;
	memset(scores, max_score, 1);
}

void SWKernel::calculate_alignment_matrix(char const * const * const read,
		char const * const * const ref, alnMat const matrix, short * const best_coordinates) {

	short best_read_pos = 0;
	short best_ref_pos = 0;
	short max_score = 0;

	short prev_row_score = 0;
	short current_row_score = 1;

	short current_row_aln = 1;

	short * scoreMat = new short[(refLength + 1) * 2]();

	for (int read_pos = 0; read_pos < readLength; ++read_pos) {

		for (int ref_pos = 0; ref_pos < refLength; ++ref_pos) {

			short up = scoreMat[prev_row_score * (refLength + 1) + ref_pos + 1];
			short diag = scoreMat[prev_row_score * (refLength + 1) + ref_pos];
			short left = scoreMat[current_row_score * (refLength + 1) + ref_pos];

			read[0][read_pos] == ref[0][ref_pos] ? diag += scoreMatch : diag += scoreMismatch;

			short cur = max(up + scoreGapRef, max(left + scoreGapRead, max(diag, 0)));

			scoreMat[current_row_score * (refLength + 1) + ref_pos + 1] = cur;

			char pointer = START;

			if (cur == 0) {
				pointer = START;
			} else if (cur == diag) {
				pointer = DIAG;
			} else if (cur == up + scoreGapRef) {
				pointer = UP;
			} else if (cur == left + scoreGapRead) {
				pointer = LEFT;
			}

			matrix[current_row_aln * (refLength + 1) + ref_pos + 1] = pointer;

			if (cur > max_score) {
				best_read_pos = read_pos;
				best_ref_pos = ref_pos;
				max_score = cur;
			}

			std::cout << cur << " ";
		}
		std::cout << std::endl;

		prev_row_score = current_row_score;
		(++current_row_score) &= 1;

		++current_row_aln;
	}

	best_coordinates[0] = best_read_pos;
	best_coordinates[1] = best_ref_pos;

}

void SWKernel::calc_alignment(char const * const * const read,
		char const * const * const ref, Alignment * const alignment) {

	alnMat matrix = new char [alnLength];
	memset(matrix, START, alnLength * sizeof(char));

	short * best_coordinates = new short[2];

	std::cout << "Score matrix" << std::endl;

	calculate_alignment_matrix(read, ref, matrix, best_coordinates);

	std::cout << "Matrix:" << std::endl;
	for (int i = 0; i < readLength + 1; ++i) {
		for (int j = 0; j < refLength + 1; ++j) {
			std::cout << matrix[i * (refLength + 1) + j] << " ";
		}
		std::cout << std::endl;
	}

	char * alignments = new char[alnLength * 2];

	int read_pos = best_coordinates[0];
	int ref_pos = best_coordinates[1];

	std::cout << "Best read: " << read_pos << std::endl << "Best ref: " << ref_pos << std::endl;

	int aln_pos = alnLength - 2;

	char backtrack = matrix[(read_pos + 1) * (refLength + 1) + ref_pos + 1];

	std::cout << "Backtrack start:\t";

	while(backtrack != START) {

		std::cout << backtrack;

		if (backtrack == UP) {
			alignments[alnLength + aln_pos] = '-';
			alignments[aln_pos] = read[0][read_pos--];
		}
		if (backtrack == LEFT) {
			alignments[aln_pos] = '-';
			alignments[alnLength + aln_pos] = ref[0][ref_pos--];
		}
		if (backtrack == DIAG) {
			alignments[aln_pos] = read[0][read_pos--];
			alignments[alnLength + aln_pos] = ref[0][ref_pos--];
		}
		backtrack = matrix[(read_pos + 1) * (refLength + 1) + ref_pos + 1];
		--aln_pos;
	}

	std::cout << "\tBacktrack end" << std::endl;

	std::cout << "==================" << std::endl;
	for (int i = 0; i < alnLength; ++i) {
		std::cout << alignments[i];
	}
	std::cout << std::endl;
	for (int i = alnLength; i < alnLength * 2; ++i) {
		std::cout << alignments[i];
	}
	std::cout << std::endl << "==================" << std::endl;
//
//
//
//	//while()
//
//	delete [] matrix; matrix = 0;
}
