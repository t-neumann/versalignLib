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
