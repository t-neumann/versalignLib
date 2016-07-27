/*
 * versalignUtil.h
 *
 *  Created on: May 25, 2016
 *      Author: tobias.neumann
 */

#ifndef VERSALIGNUTIL_H_
#define VERSALIGNUTIL_H_

#include <string.h>
#include <sstream>
#include <immintrin.h>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <iostream>

#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

#ifdef WIN32   // Windows system specific
#include <windows.h>
#else          // Unix based system specific
#include <sys/time.h>
#endif

using std::ifstream;
using std::cout;
using std::getline;
using std::string;
using std::vector;

#include "AlignmentParameters.h"
#include "AlignmentLogger.h"

/*############
 * CLASSES
 *###########*/

class FastaProvider {

public:
	FastaProvider() {
	}

	vector<const char *> parse_fasta(string const & filename) {

		vector<const char*> values;

		ifstream in(filename.c_str());

		if (!in.good()) {
			std::cerr << "Error opening " + filename << ".\n";
		} else {
			std::string line, name, content;
			while (getline(in, line).good()) {
				if (line.empty() || line[0] == '>') {
					if (!name.empty()) {
						//char * charSeq = new char[content.length() + 1];
						//strncpy(charSeq, content.c_str(), content.length() + 1);
						//values.push_back(charSeq);
						values.push_back(strdup(content.c_str()));
						name.clear();
					}
					if (!line.empty()) {
						name = line.substr(1);
					}
					content.clear();
				} else if (!name.empty()) {
					if (line.find(' ') != std::string::npos) {
						name.clear();
						content.clear();
					} else {
						content += line;
					}
				}
			}
			if (!name.empty()) {
				//char * charSeq = new char[content.length() + 1];
				//strncpy(charSeq, content.c_str(), content.length() + 1);
				//values.push_back(charSeq);
				values.push_back(strdup(content.c_str()));
			}
		}
		return values;
	}



	~FastaProvider() {
	}
};

//////////////////////////////////////////////////////////////////////////////
// Timer.h
// =======
// High Resolution Timer.
// This timer is able to measure the elapsed time with 1 micro-second accuracy
// in both Windows, Linux and Unix system
//
//  AUTHOR: Song Ho Ahn (song.ahn@gmail.com)
// CREATED: 2003-01-13
// UPDATED: 2006-01-13
//
// Copyright (c) 2003 Song Ho Ahn
//////////////////////////////////////////////////////////////////////////////

class Timer
{
public:
    Timer();                                    // default constructor
    ~Timer();                                   // default destructor

    void   start();                             // start timer
    void   stop();                              // stop the timer
    double getElapsedTime();                    // get elapsed time in second
    double getElapsedTimeInSec();               // get elapsed time in second (same as getElapsedTime)
    double getElapsedTimeInMilliSec();          // get elapsed time in milli-second
    double getElapsedTimeInMicroSec();          // get elapsed time in micro-second


protected:


private:
    double startTimeInMicroSec;                 // starting time in micro-second
    double endTimeInMicroSec;                   // ending time in micro-second
    int    stopped;                             // stop flag
#ifdef WIN32
    LARGE_INTEGER frequency;                    // ticks per second
    LARGE_INTEGER startCount;                   //
    LARGE_INTEGER endCount;                     //
#else
    timeval startCount;                         //
    timeval endCount;                           //
#endif
};


/*############
 * FUNCTIONS
 *###########*/

void * DLL_function_retreival(int const dll, char const * const name, bool required = true);

int const DLL_init(char const * const filename, AlignmentParameters * parameters, AlignmentLogger * logger);

size_t pad(char const * * strings, int const & n, char const & pad);

bool check_avx2_support ();

/*############
 * TEMPLATES
 *###########*/

template <typename T>
std::string __m128i_toString(const __m128i var) {
    std::stringstream sstr;
    const T* values = (const T*) &var;
    if (sizeof(T) == 1) {
        for (unsigned int i = 0; i < sizeof(__m128i); i++) {
            sstr << (int) values[i] << " ";
        }
    } else {
        for (unsigned int i = 0; i < sizeof(__m128i) / sizeof(T); i++) {
            sstr << values[i] << " ";
        }
    }
    return sstr.str();
}

#endif /* VERSALIGNUTIL_H_ */
