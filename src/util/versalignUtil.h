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
#include <emmintrin.h>

size_t pad(char const * * strings, int const & n, char const & pad);

template <typename T> std::string __m128i_toString(const __m128i var);

#endif /* VERSALIGNUTIL_H_ */
