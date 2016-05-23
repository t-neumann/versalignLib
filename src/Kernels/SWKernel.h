/*
 * SWKernel.h
 *
 *  Created on: May 23, 2016
 *      Author: tobias.neumann
 */

#ifndef SWKERNEL_H_
#define SWKERNEL_H_

#include "AlignmentKernel.h"

#include <cstring>
#include <cmath>
#include <iostream>
#include <xmmintrin.h>
#include <stdlib.h>
#include <cstring>
#include <ctime>

class SWKernel: public AlignmentKernel {

public:

	SWKernel();
	virtual ~SWKernel();

	virtual float score_alignment(char const * read, float const & read_length,
			char const * reference, float const & reference_length,
			float const & gap_read, float const & gap_ref, float const & match,
			float const & mismatch);
	virtual float score_alignment_corridor(char const * read, float const & read_length,
			char const * reference, float const & reference_length,
			float const & gap_read, float const & gap_ref, float const & match,
			float const & mismatch, float const & corridor_width);
};

#endif /* SWKERNEL_H_ */
