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

			diag += base_score[char_to_score[read[0][read_pos]]][char_to_score[ref[0][ref_pos]]];

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

void SWKernel::score_alignment_needleman_wunsch(char const * const * const read,
		char const * const * const ref, short * const scores) {

	short * matrix = new short[(refLength + 1) * 2]();

	int prev_row = 0;
	int current_row = 1;

	short globalMax = 0;

	for (int read_pos = 0; read_pos < readLength; ++read_pos) {

		for (int ref_pos = 0; ref_pos < refLength; ++ref_pos) {

			short up = matrix[prev_row * (refLength + 1) + ref_pos + 1];
			short diag = matrix[prev_row * (refLength + 1) + ref_pos];
			short left = matrix[current_row * (refLength + 1) + ref_pos];

			diag += base_score[char_to_score[read[0][read_pos]]][char_to_score[ref[0][ref_pos]]];

			short cur = max(up + scoreGapRef, max(left + scoreGapRead, diag));

			matrix[current_row * (refLength + 1) + ref_pos + 1] = cur;

			//std::cout << cur << " ";
		}
		//std::cout << std::endl;

		globalMax = max(globalMax, matrix[current_row * (refLength + 1) + refLength]);

		prev_row = current_row;
		(++current_row) &= 1;
	}

	for (int ref_pos = 0; ref_pos < refLength + 1; ++ref_pos) {
		globalMax = max(globalMax, matrix[prev_row * (refLength + 1) + ref_pos]);
	}

	memset(scores, globalMax, 1);

	delete [] matrix; matrix = 0;
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

			diag += base_score[char_to_score[read[0][read_pos]]][char_to_score[ref[0][ref_pos]]];

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

	delete []scoreMat; scoreMat = 0;

	best_coordinates[0] = best_read_pos;
	best_coordinates[1] = best_ref_pos;
}

void SWKernel::calculate_alignment_matrix_needleman_wunsch(char const * const * const read,
		char const * const * const ref, alnMat const matrix, short * const best_coordinates) {

	//for (int ref_pos = 0; ref_pos < refLength; ++ref_pos) {
	//	matrix[ref_pos + 1] = LEFT;
	//}

	//short max_read_pos = SHRT_MIN;
	//short max_ref_pos = SHRT_MIN;

	short max_read_pos = readLength - 1;
	short max_ref_pos = refLength - 1;

	short prev_row_score = 0;
	short current_row_score = 1;

	short current_row_aln = 1;

	short * scoreMat = new short[(refLength + 1) * 2]();

	short rowMax = SHRT_MIN;
	//short colMax = SHRT_MIN;

//	short globalRowMax = SHRT_MIN;
	short globalRowMaxIndex = -1;

	short rowMaxIndex = 0;
	//short colMaxIndex = 0;

	for (int read_pos = 0; read_pos < readLength + 1; ++read_pos) {
//		std::cout << scoreMat[read_pos] << " ";
	}
//	std::cout << std::endl;

	for (int read_pos = 0; read_pos < readLength; ++read_pos) {

		matrix[current_row_aln * (refLength + 1)] = UP;
		scoreMat[current_row_score * (refLength + 1)] = (read_pos + 1) * scoreGapRef;

		// 0 score means forbidden char
		if (max_read_pos == readLength - 1 && char_to_score[read[0][read_pos]] == 0) {
			max_read_pos = read_pos - 1;
		}

		// Save previous row max if read ends prematurely
		if (max_read_pos + 1 == read_pos) {
//			globalRowMax = rowMax;
			globalRowMaxIndex = rowMaxIndex;
		}

		rowMax = scoreMat[current_row_score * (refLength + 1)];
		rowMaxIndex = 0;

//		std::cout << (read_pos + 1) * scoreGapRef << " ";

		for (int ref_pos = 0; ref_pos < refLength; ++ref_pos) {

			short up = scoreMat[prev_row_score * (refLength + 1) + ref_pos + 1];
			short diag = scoreMat[prev_row_score * (refLength + 1) + ref_pos];
			short left = scoreMat[current_row_score * (refLength + 1) + ref_pos];

			diag += base_score[char_to_score[read[0][read_pos]]][char_to_score[ref[0][ref_pos]]];

			short cur = max(up + scoreGapRef, max(left + scoreGapRead, diag));

			scoreMat[current_row_score * (refLength + 1) + ref_pos + 1] = cur;

			char pointer = START;

			if (cur == diag) {
				pointer = DIAG;
			} else if (cur == up + scoreGapRef) {
				pointer = UP;
			} else if (cur == left + scoreGapRead) {
				pointer = LEFT;
			}

			if (max_ref_pos == refLength - 1 && char_to_score[ref[0][ref_pos]] == 0) {
				//std::cout << std::endl << "Max ref is " << ref_pos - 1 << std::endl;
				max_ref_pos = ref_pos - 1;
			}

			if (cur > rowMax) {
				//std::cout << std::endl << "Max row is " << cur << " Max pos is " << ref_pos << std::endl;
				rowMax = cur;
				rowMaxIndex = ref_pos;
			}

			matrix[current_row_aln * (refLength + 1) + ref_pos + 1] = pointer;

//			std::cout << cur << " ";
		}
//		std::cout << std::endl;

		prev_row_score = current_row_score;
		(++current_row_score) &= 1;

		++current_row_aln;
	}

	delete []scoreMat; scoreMat = 0;

	best_coordinates[0] = max_read_pos;

//	std::cout << "Max_ref_pos: " << max_ref_pos << std::endl << "globalRowMaxIndex: " << globalRowMaxIndex << std::endl << "rowMaxIndex: " << rowMaxIndex << std::endl;


	if (globalRowMaxIndex < 0 ) {
		//globalRowMax = rowMax;
		globalRowMaxIndex = rowMaxIndex;
	}

//	if (max_ref_pos == refLength - 1) {
//		best_coordinates[1] = globalRowMaxIndex;
//	} else {
//		best_coordinates[1] = max_ref_pos;
//	}
	best_coordinates[1] = std::min(max_ref_pos, globalRowMaxIndex);

}

void SWKernel::calc_alignment(char const * const * const read,
		char const * const * const ref, Alignment * const alignment) {

	std::cout << "Aln Length:\t" << alnLength << std::endl;

	alnMat matrix = new char [(refLength + 1) * (readLength + 1)];
	memset(matrix, START, (refLength + 1) * (readLength + 1) * sizeof(char));

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

	alignment->read = new char[alnLength];
	alignment->ref = new char[alnLength];

	memcpy(alignment->read, alignments, alnLength * sizeof(char));
	memcpy(alignment->ref, alignments + alnLength, alnLength * sizeof(char));

	alignment->readStart = aln_pos + 1;
	alignment->refStart = aln_pos + 1;

	alignment->readEnd = alnLength - 1;
	alignment->refEnd = alnLength - 1;

	std::cout << "\tBacktrack end" << std::endl;

	delete [] alignments; alignments = 0;
	delete [] best_coordinates; best_coordinates = 0;
	delete [] matrix; matrix = 0;
}

void SWKernel::calc_alignment_needleman_wunsch(char const * const * const read,
		char const * const * const ref, Alignment * const alignment) {

	std::cout << "Aln Length:\t" << alnLength << std::endl;

	alnMat matrix = new char [(refLength + 1) * (readLength + 1)];
	memset(matrix, START, (refLength + 1) * (readLength + 1) * sizeof(char));

	short * max_coordinates = new short[2];

	std::cout << "Score matrix" << std::endl;

	calculate_alignment_matrix_needleman_wunsch(read, ref, matrix, max_coordinates);

	std::cout << "Matrix:" << std::endl;
	for (int i = 0; i < readLength + 1; ++i) {
		for (int j = 0; j < refLength + 1; ++j) {
//			std::cout << matrix[i * (refLength + 1) + j] << " ";
		}
//		std::cout << std::endl;
	}

	char * alignments = new char[alnLength * 2];

	int read_pos = max_coordinates[0];
	int ref_pos = max_coordinates[1];

	std::cout << "Best read: " << read_pos << std::endl << "Best ref: " << ref_pos << std::endl;

	int aln_pos = alnLength - 2;

	char backtrack = matrix[(read_pos + 1) * (refLength + 1) + ref_pos + 1];

	std::cout << "Backtrack start:\t";

	alignments[alnLength - 1] = '\0';
	alignments[2 * alnLength - 1] = '\0';

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

	alignment->read = new char[alnLength];
	alignment->ref = new char[alnLength];

	memcpy(alignment->read, alignments, alnLength * sizeof(char));
	memcpy(alignment->ref, alignments + alnLength, alnLength * sizeof(char));

	alignment->readStart = aln_pos + 1;
	alignment->refStart = aln_pos + 1;

	alignment->readEnd = alnLength - 1;
	alignment->refEnd = alnLength - 1;

	std::cout << "\tBacktrack end" << std::endl;

	delete [] alignments; alignments = 0;
	delete [] max_coordinates; max_coordinates = 0;
	delete [] matrix; matrix = 0;
}
