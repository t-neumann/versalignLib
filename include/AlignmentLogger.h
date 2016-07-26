/*
 * AlignmentLogger.h
 *
 *  Created on: Jul 26, 2016
 *      Author: tobias.neumann
 */

#ifndef ALIGNMENTLOGGER_H
#define ALIGNMENTLOGGER_H

#include <stdio.h>

class AlignmentLogger {
	public:
		virtual void log(int const level, char const * const main, char const * const msg, size_t const & arg_num = 0, ...) = 0;
		virtual ~AlignmentLogger() {}
};

typedef void (*fp_set_logger)(AlignmentLogger const *);

extern AlignmentLogger* _logger;
#define Logger (*_logger)

#endif /* ALIGNMENTLOGGER_H */
