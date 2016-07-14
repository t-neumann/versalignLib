/*
 * AlignmentParameters.h
 *
 *  Created on: Jul 13, 2016
 *      Author: tobias.neumann
 */

#ifndef INCLUDE_ALIGNMENTPARAMETERS_H
#define INCLUDE_ALIGNMENTPARAMETERS_H

class AlignmentParameters {
public:
	virtual int param_int(char const * const key) = 0;
	virtual bool has_key(char const * const key) = 0;

	virtual ~AlignmentParameters() {};
};

typedef void (*fp_set_parameters)(AlignmentParameters const *);

extern AlignmentParameters* _parameters;
#define Parameters (*_parameters)

#endif /* INCLUDE_ALIGNMENTPARAMETERS_H */
