/*
 * dll_handling.h
 *
 *  Created on: Jul 14, 2016
 *      Author: tobias.neumann
 */

#ifndef DLL_HANDLING_H
#define DLL_HANDLING_H

#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

#include "AlignmentParameters.h"

void * DLL_function_retreival(int const dll, char const * const name, bool required = true);

int const DLL_init(char const * const filename, AlignmentParameters * parameters);

#endif /* DLL_HANDLING_H */
