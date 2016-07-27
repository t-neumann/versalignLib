cmake_minimum_required (VERSION 2.8.11)

if (AVX_FLAG AND AVX2_FLAG)
	message("Compiling AVX Kernel")
	add_library(AVXKernel SHARED AVXKernel.cpp AVXKernel_dllexport.cpp)	
else()
	message("No AVX support, ignoring AVX Kernel")
endif()
    
add_library(SSEKernel SHARED SSEKernel.cpp SSEKernel_dllexport.cpp)
