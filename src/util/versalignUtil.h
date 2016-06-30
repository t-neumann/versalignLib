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

size_t pad(char const * * strings, int const & n, char const & pad);

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
