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
#include <sstream>
#include <string>
#include <stdarg.h>

class CustomLogger: public AlignmentLogger {
public:
	virtual void log(int const level, char const * const main,
			char const * const msg, size_t const & arg_num = 0, ...) {

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

		sstream << severity << "\t[" << main << "]\t" << msg;

		if (arg_num > 0) {
			va_list args;
			va_start(args, arg_num);
			for (int i = 0; i < arg_num; ++i) {
				sstream << std::endl << severity << "\t[" << main << "]\t" << va_arg(args, const char *);
			}
		}

		sstream << std::endl;

		std::cerr << sstream.str();

		clear_buffer();
	}
private:
	std::stringstream sstream;

	void clear_buffer() {
		sstream.clear();
		sstream.str(std::string());
	}
};

#endif /* CUSTOMLOGGER_H */
