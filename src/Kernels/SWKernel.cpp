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

xy_coordinates SWKernel::calculate_alignment_matrix(char const * const * const read,
		char const * const * const ref, mat_element * const matrix) {

	short best_read_pos = 0;
	short best_ref_pos = 0;
	short max_score = 0;

	int prev_row = 0;
	int current_row = 1;

	for (int read_pos = 0; read_pos < readLength; ++read_pos) {

		for (int ref_pos = 0; ref_pos < refLength; ++ref_pos) {

			mat_element up = matrix[prev_row * (refLength + 1) + ref_pos + 1];
			mat_element diag = matrix[prev_row * (refLength + 1) + ref_pos];
			mat_element left = matrix[current_row * (refLength + 1) + ref_pos];

			short compareScore = 0;
			read[0][read_pos] == ref[0][ref_pos] ? compareScore = diag.score + scoreMatch : compareScore = diag.score + scoreMismatch;

			short curScore = max(up.score + scoreGapRef, max(left.score + scoreGapRead, max(compareScore, 0)));

			matrix[current_row * (refLength + 1) + ref_pos + 1].score = curScore;

			if(up.score + scoreGapRef == curScore) {
				matrix[current_row * (refLength + 1) + ref_pos].path = upwards;
			} else if (left.score + scoreGapRead == curScore) {
				matrix[current_row * (refLength + 1) + ref_pos].path = leftwards;
			} else if (diag.score + compareScore == curScore) {
				matrix[current_row * (refLength + 1) + ref_pos].path = diagonal;
			}

			if (curScore > max_score) {
				best_read_pos = read_pos;
				best_ref_pos = ref_pos;
				max_score = curScore;
			}
		}
		prev_row = current_row;
		++current_row;
	}

	xy_coordinates best_score;
	best_score.x = best_read_pos;
	best_score.y = best_ref_pos;

	return best_score;
}

void SWKernel::calc_alignment(char const * const * const read,
		char const * const * const ref, char * const * const aligned_read,
		char * const * const aligned_ref) {

	mat_element * matrix = new mat_element[(refLength + 1) * (readLength + 1)]();

	xy_coordinates best_score = calculate_alignment_matrix(read, ref, matrix);

	int read_pos = best_score.x;
	int ref_pos = best_score.y;

	mat_element backtrack = matrix[(read_pos + 1) * (refLength + 1) + ref_pos + 1];

	while(backtrack.path != start) {
		std::cout << backtrack.path;

		if (backtrack.path == upwards) {
			--read_pos;
		}
		if (backtrack.path == leftwards) {
			--ref_pos;
		}
		if (backtrack.path == diagonal) {
			--ref_pos;
			--read_pos;
		}
		backtrack = matrix[(read_pos + 1) * (refLength + 1) + ref_pos + 1];
	}

	std::cout << "Best scores:\t" << best_score.x << "\t" << best_score.y << "\t" << backtrack.score << std::endl;

	//while()

	delete [] matrix; matrix = 0;
}
