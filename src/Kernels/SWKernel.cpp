/*
 * SWAligner.cpp
 *
 *  Created on: May 23, 2016
 *      Author: tobias.neumann
 */

#include "SWKernel.h"

using std::cout;
using std::endl;

// Helper declaration

float get_score(char const & a, char const & b, float const & match,
		float const & mismatch);

void print_row(char const & a, float * row, float const & reference_length);

// Class method implementation

SWKernel::SWKernel() {
	// TODO Auto-generated constructor stub

}

SWKernel::~SWKernel() {
	// TODO Auto-generated destructor stub
}

// Main functions

float SWKernel::score_alignment(char const * read, float const & read_length,
		char const * reference, float const & reference_length,
		float const & gap_read, float const & gap_ref, float const & match,
		float const & mismatch) {

	cout << "  - ";
	for (int i = 0; i < reference_length; ++i) {
		cout << reference[i] << " ";
	}
	cout << endl;

	float max_score = 0.0f;

	float * pPrev = new float[(int) reference_length + 1];
	float * pNext = new float[(int) reference_length + 1];

	for (int i = 0; i <= (int) reference_length; ++i) {
		pPrev[i] = 0.0f;
	}

	print_row('-', pPrev, reference_length);

	for (int read_index = 0; read_index < (int) read_length; ++read_index) {
		for (int reference_index = 0; reference_index <= (int) reference_length;
				++reference_index) {
			pNext[reference_index] = std::max(pPrev[reference_index] + gap_ref,
					0.0f);
			if (reference_index > 0) {
				pNext[reference_index] = std::max(pNext[reference_index],
						pNext[reference_index - 1] + gap_read);
				pNext[reference_index] = std::max(pNext[reference_index],
						pPrev[reference_index - 1]
								+ get_score(read[read_index],
										reference[reference_index - 1], match,
										mismatch));
			}

			if (pNext[reference_index] > max_score) {
				max_score = pNext[reference_index];
			}
		}
		print_row(read[read_index], pNext, reference_length);
		delete[] pPrev;
		pPrev = pNext;
		pNext = new float[(int) reference_length + 1];
	}
	delete []pNext; pNext = 0;

	return max_score;
}

float SWKernel::score_alignment_corridor(char const * read,
		float const & read_length, char const * reference,
		float const & reference_length, float const & gap_read,
		float const & gap_ref, float const & match, float const & mismatch,
		float const & corridor_width) {

	float max_score = 0.0f;

	float * pPrev = new float[(int) corridor_width];
	float * pNext = new float[(int) corridor_width];

	for (int i = 0; i < (int) corridor_width; ++i) {
		pPrev[i] = 0.0f;
	}

	for (int read_index = 0; read_index < (int) read_length; ++read_index) {
		for (int reference_index = 0; reference_index < (float) corridor_width;
				++reference_index) {

			pNext[reference_index] = std::max(0.0f,
					pPrev[reference_index]
							+ get_score(read[read_index],
									reference[reference_index - 1 + read_index],
									match, mismatch));

			if (reference_index < corridor_width - 1) {
				pNext[reference_index] = std::max(
						pPrev[reference_index + 1] + gap_ref,
						pNext[reference_index]);
			}

			if (reference_index > 0) {
				pNext[reference_index] = std::max(pNext[reference_index],
						pNext[reference_index - 1] + gap_read);
			}

			if (pNext[reference_index] > max_score) {
				max_score = pNext[reference_index];
			}
		}

		delete[] pPrev;
		pPrev = pNext;
		pNext = new float[(int) corridor_width];
	}

	return max_score;
}

float get_score(char const & a, char const & b, float const & match,
		float const & mismatch) {
	return a == b ? match : mismatch;
}

void print_row(char const & a, float * row, float const & reference_length) {

	for (int i = 0; i <= reference_length; ++i) {
		std::cout << row[i] << " ";
	}
	std::cout << std::endl;
}

