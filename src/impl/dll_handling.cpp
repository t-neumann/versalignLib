#include "dll_handling.h"

#include <map>
#include <iostream>

std::map<int, void *> DLLmap;

void * DLL_function_retreival(int const dll, char const * const name,
		bool required) {
	void * sdll = DLLmap[dll];
	void * pf = dlsym(sdll, name);
	if (required && (pf == 0)) {
		std::cerr << "COULD NOT FIND FUNCTION " << name << std::endl;
	}
	return pf;
}

int const DLL_init(char const * const filename,
		AlignmentParameters * parameters, AlignmentLogger * logger) {
	try {
		void * sdl_lib = dlopen(filename, RTLD_LAZY);

		if (!sdl_lib) {
			std::cerr << "Failed loading DLL: " << filename << std::endl;
			return -1;
		}

		int index = DLLmap.size();
		DLLmap[index] = sdl_lib;

		fp_set_parameters set_parameters =
				(fp_set_parameters) DLL_function_retreival(index,
						"set_parameters");
		if (set_parameters != 0) {
			set_parameters(parameters);
		}

		fp_set_logger set_logger = (fp_set_logger) DLL_function_retreival(index,
				"set_logger");
		if (set_logger != 0) {
			set_logger(logger);
		}

		return index;
	} catch (int err) {
		std::cerr << "Error opening DLL: " << filename << std::endl;
	}
	return 0;
}
