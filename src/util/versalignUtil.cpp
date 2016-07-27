/*

 * versalignUtil.cpp
 *
 *  Created on: May 25, 2016
 *      Author: tobias.neumann
 */

#include "versalignUtil.h"

std::map<int, void *> DLLmap;

inline int max(int a, int b) {
  return a > b ? a : b;
}

size_t pad(char const * * strings, int const & n, char const & pad) {
	size_t max_length = 0;
	for (int i = 0; i < n; ++i) {
		max_length = max(max_length, strlen(strings[i]));
	}
	// Add one for null character termination
	//++max_length;
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

/*
 * 	Source downloaded on July 15, 2016:
 * 	https://software.intel.com/en-us/articles/how-to-detect-new-instruction-support-in-the-4th-generation-intel-core-processor-family
 *
 *	Modified by Tobias Neumann
 *	Email: tobias.neumann.at@gmail.com
 */

// Switch to __builtin_cpu_init()
// and then  __builtin_cpu_supports("avx2") for future GCC/clang versions (>= 4.8.5)
// https://gcc.gnu.org/onlinedocs/gcc-4.8.5/gcc/X86-Built-in-Functions.html


#if defined(__INTEL_COMPILER) && (__INTEL_COMPILER >= 1300)

#include <immintrin.h>

int check_4th_gen_intel_core_features()
{
    const int the_4th_gen_features =
        (_FEATURE_AVX2 | _FEATURE_FMA | _FEATURE_BMI | _FEATURE_LZCNT | _FEATURE_MOVBE);
    return _may_i_use_cpu_feature( the_4th_gen_features );
}

#else /* non-Intel compiler */

#include <stdint.h>
#if defined(_MSC_VER)
# include <intrin.h>
#endif

void run_cpuid(uint32_t eax, uint32_t ecx, uint32_t* abcd)
{
#if defined(_MSC_VER)
    __cpuidex(abcd, eax, ecx);
#else
    uint32_t ebx, edx;
# if defined( __i386__ ) && defined ( __PIC__ )
     /* in case of PIC under 32-bit EBX cannot be clobbered */
    __asm__ ( "movl %%ebx, %%edi \n\t cpuid \n\t xchgl %%ebx, %%edi" : "=D" (ebx),
# else
    __asm__ ( "cpuid" : "+b" (ebx),
# endif
              "+a" (eax), "+c" (ecx), "=d" (edx) );
    abcd[0] = eax; abcd[1] = ebx; abcd[2] = ecx; abcd[3] = edx;
#endif
}

int check_xcr0_ymm()
{
    uint32_t xcr0;
#if defined(_MSC_VER)
    xcr0 = (uint32_t)_xgetbv(0);  /* min VS2010 SP1 compiler is required */
#else
    __asm__ ("xgetbv" : "=a" (xcr0) : "c" (0) : "%edx" );
#endif
    return ((xcr0 & 6) == 6); /* checking if xmm and ymm state are enabled in XCR0 */
}

int check_4th_gen_intel_core_features() {
    uint32_t abcd[4];
    uint32_t fma_movbe_osxsave_mask = ((1 << 12) | (1 << 22) | (1 << 27));
    uint32_t avx2_bmi12_mask = (1 << 5) | (1 << 3) | (1 << 8);

    /* CPUID.(EAX=01H, ECX=0H):ECX.FMA[bit 12]==1   &&
       CPUID.(EAX=01H, ECX=0H):ECX.MOVBE[bit 22]==1 &&
       CPUID.(EAX=01H, ECX=0H):ECX.OSXSAVE[bit 27]==1 */
    run_cpuid( 1, 0, abcd );
    if ( (abcd[2] & fma_movbe_osxsave_mask) != fma_movbe_osxsave_mask )
        return 0;

    if ( ! check_xcr0_ymm() )
        return 0;

    /*  CPUID.(EAX=07H, ECX=0H):EBX.AVX2[bit 5]==1  &&
        CPUID.(EAX=07H, ECX=0H):EBX.BMI1[bit 3]==1  &&
        CPUID.(EAX=07H, ECX=0H):EBX.BMI2[bit 8]==1  */
    run_cpuid( 7, 0, abcd );
    if ( (abcd[1] & avx2_bmi12_mask) != avx2_bmi12_mask )
        return 0;

    /* CPUID.(EAX=80000001H):ECX.LZCNT[bit 5]==1 */
    run_cpuid( 0x80000001, 0, abcd );
    if ( (abcd[2] & (1 << 5)) == 0)
        return 0;

    return 1;
}

#endif /* non-Intel compiler */

static int checked = -1;
static bool the_4th_gen_features_available = false;

bool check_avx2_support() {
    /* test is performed once */
    if (checked < 0 ) {
    	checked = check_4th_gen_intel_core_features();
    	if (checked > 0) {
    		the_4th_gen_features_available = true;
    	}
    }
    return the_4th_gen_features_available;
}
