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
	char * read = 0;
	char * ref = 0;
	short readStart;
	short readEnd;
	short refStart;
	short refEnd;

	~Alignment() {
		if (read != 0) delete[] read;
		if (ref != 0) delete[] ref;
	}
};

/*
	Options:

	int alignment_algorithm = opt & 0xF;
			-> 0 .. Smith-Waterman
			-> 1 .. Needleman-Wunsch
*/

class AlignmentKernel {

public:

	virtual ~AlignmentKernel() {}

	virtual void score_alignments(int const & opt, int const & aln_number, char const * const * const reads,
			char const * const * const refs, short * const scores) = 0;
	virtual void compute_alignments(int const & opt, int const & aln_number, char const * const * const reads,
			char const * const * const refs, Alignment * const alignments) = 0;
};

typedef AlignmentKernel * (*fp_load_alignment_kernel)();
typedef void (*fp_delete_alignment_kernel)(AlignmentKernel*);



#endif /* ALIGNMENTKERNEL_H */
