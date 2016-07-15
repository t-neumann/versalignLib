MESSAGE( STATUS "CMAKE_CURRENT_BINARY_DIR: " ${CMAKE_CURRENT_BINARY_DIR} )
MESSAGE( STATUS "CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR} )

if (AVX_FLAG AND AVX2_FLAG)
	message("Compiling AVX Kernel")
	add_library(AVXKernel SHARED AVXKernel.cpp AVXKernel_dllexport.cpp)	
else()
	message("No AVX support, ignoring AVX Kernel")
endif()
    
add_library(SSEKernel SHARED SSEKernel.cpp SSEKernel_dllexport.cpp)
