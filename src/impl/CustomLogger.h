/*
 * CustomLogger.h
 *
 *  Created on: Jul 26, 2016
 *      Author: tobias.neumann
 */

#ifndef CUSTOMLOGGER_H
#define CUSTOMLOGGER_H

#include "AlignmentLogger.h"
#include <iostream>

class CustomLogger : public AlignmentLogger {
	public:
		void log(int const level, char const * const head, char const * const msg) const {

			char const * severity = "ERROR";
			switch (level) {
				case 0:
					severity = "INFO";
					break;
				case 1:
					severity = "WARNING";
					break;
				case 3:
					severity = "DRASTIC";
					break;
				default:
					break;
			}

			std::cerr << "[" << severity << "]\t" << head << ":\t" << msg << std::endl;
		}
};

#endif /* CUSTOMLOGGER_H */
