/*
 * SSEKernel.cpp
 *
 *  Created on: May 23, 2016
 *      Author: tobias.neumann
 */

#include "SSEKernel.h"

SSEKernel::SSEKernel() {
	// TODO Auto-generated constructor stub

}

SSEKernel::~SSEKernel() {
	// TODO Auto-generated destructor stub
}

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
