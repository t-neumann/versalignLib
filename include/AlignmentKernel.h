/*
 * AlignmentKernel.h
 *
 *  Created on: May 23, 2016
 *      Author: Tobias Neumann
 *      Email: tobias.neumann.at@gmail.com
 */

#ifndef ALIGNMENTKERNEL_H
#define ALIGNMENTKERNEL_H

struct Alignment {
	char * read;
	char * ref;
	short readStart;
	short readEnd;
	short refStart;
	short refEnd;
};

class AlignmentKernel {

public:

	virtual ~AlignmentKernel() {}

	virtual void score_alignment(int const & opt, int const & aln_number, char const * const * const reads,
			char const * const * const refs, short * const scores) = 0;
	virtual void compute_alignment(int const & opt, int const & aln_number, char const * const * const reads,
			char const * const * const refs, Alignment * const alignments) = 0;
};

typedef AlignmentKernel * (*load_alignment_kernel)();
typedef void (*delete_alignment_kernel)(AlignmentKernel*);



#endif /* ALIGNMENTKERNEL_H */
