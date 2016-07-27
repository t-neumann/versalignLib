MESSAGE( STATUS "CMAKE_CURRENT_BINARY_DIR: " ${CMAKE_CURRENT_BINARY_DIR} )
MESSAGE( STATUS "CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR} )

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/src/cmake")

include_directories (.)
include_directories ("${PROJECT_SOURCE_DIR}/include")

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

include(FindOpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_subdirectory (util)
add_subdirectory (Kernels)
add_subdirectory (impl)