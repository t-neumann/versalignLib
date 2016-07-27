MESSAGE( STATUS "CMAKE_CURRENT_BINARY_DIR: " ${CMAKE_CURRENT_BINARY_DIR} )
MESSAGE( STATUS "CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR} )

# -O3 apparently doesn't work with SSE instructions ("Illegal instruction") on gcc/g++

SET(CMAKE_C_FLAGS_RELEASE "-O2 -DNDEBUG")
SET(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")

if (AVX_FLAG AND AVX2_FLAG)
	message("Compiling AVX Kernel")
	add_library(AVXKernel SHARED AVXKernel.cpp AVXKernel_dllexport.cpp)	
else()
	message("No AVX support, ignoring AVX Kernel")
endif()
    
add_library(SSEKernel SHARED SSEKernel.cpp SSEKernel_dllexport.cpp)
