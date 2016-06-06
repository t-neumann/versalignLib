/*

 * versalignUtil.cpp
 *
 *  Created on: May 25, 2016
 *      Author: tobias.neumann
 */

#include "versalignUtil.h"

inline int max(int a, int b) {
  return a > b ? a : b;
}

size_t pad(char const * * strings, int const & n, char const & pad) {
	size_t max_length = 0;
	for (int i = 0; i < n; ++i) {
		max_length = max(max_length, strlen(strings[i]));
	}
	for (int i = 0; i < n; ++i) {
		char * tmp = new char[max_length];
		size_t length = strlen(strings[i]);
		for (size_t j = 0; j < max_length; ++j) {
			tmp[j] = (j < length) ? strings[i][j] : pad;
		}
		strings[i] = tmp;
	}
	return max_length;
}

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
