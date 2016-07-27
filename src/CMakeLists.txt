cmake_minimum_required (VERSION 2.8.11)

###################
# Setup environment
###################

# Modules

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/src/cmake")

# Includes

include_directories (.)
include_directories ("${PROJECT_SOURCE_DIR}/include")

# Check for AVX support

include(CheckCompilerFlag)

CHECK_COMPILER_FLAG(-mavx AVX_FLAG)
CHECK_COMPILER_FLAG(-mavx2 AVX2_FLAG)

if (AVX_FLAG AND AVX2_FLAG)
	message("AVX2 support detected")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx -mavx2")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx -mavx2")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mavx -mavx2")
else()
	message("No AVX2 support detected")
endif()

# Check for OpenMP support

include(FindOpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Set to C++11 version

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# -O3 apparently doesn't work with SSE instructions ("Illegal instruction") on gcc/g++

SET(CMAKE_C_FLAGS_RELEASE "-O2 -DNDEBUG")
SET(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")

add_subdirectory (util)
add_subdirectory (Kernels)
add_subdirectory (impl)