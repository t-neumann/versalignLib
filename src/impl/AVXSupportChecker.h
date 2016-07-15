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

#ifndef CHECKAVXSUPPORT_H
#define CHECKAVXSUPPORT_H

bool check_avx2_support ();

#endif /* CHECKAVXSUPPORT_H */
