/*
 * AlignmentLogger.h
 *
 *  Created on: Jul 26, 2016
 *      Author: tobias.neumann
 */

#ifndef ALIGNMENTLOGGER_H
#define ALIGNMENTLOGGER_H

class AlignmentLogger {
	public:
		virtual void log(int const lvl, char const * const title, char const * const msg) const = 0;
		virtual ~AlignmentLogger() {}
};

typedef void (*fp_set_logger)(AlignmentLogger const *);

extern AlignmentLogger* _logger;
#define Logger (*_logger)

#endif /* ALIGNMENTLOGGER_H */
