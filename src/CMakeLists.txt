MESSAGE( STATUS "CMAKE_CURRENT_BINARY_DIR: " ${CMAKE_CURRENT_BINARY_DIR} )
MESSAGE( STATUS "CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR} )

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/src/cmake")

include_directories (.)
include_directories ("${PROJECT_SOURCE_DIR}/include")

FIND_PATH(OPENCL_INCLUDE_DIR
	NAMES
		CL/cl.h OpenCL/cl.h
	PATHS
		${PROJECT_SOURCE_DIR}/include
)

message(STATUS "This is include dir " ${OPENCL_INCLUDE_DIR})

IF(CMAKE_SIZEOF_VOID_P EQUAL 4)
	SET(OPENCL_LIB_SEARCH_PATH
		${OPENCL_LIB_SEARCH_PATH}
		${PROJECT_SOURCE_DIR}/opencl-AMD-sdk/x86/lib)
ELSEIF(CMAKE_SIZEOF_VOID_P EQUAL 8)
  SET(OPENCL_LIB_SEARCH_PATH
		${OPENCL_LIB_SEARCH_PATH}
		${PROJECT_SOURCE_DIR}/opencl-AMD-sdk/x86_64/lib)
ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 4)


FIND_LIBRARY(
    OPENCL_LIBRARY
    NAMES OpenCL
    PATHS ${OPENCL_LIB_SEARCH_PATH})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  OpenCL
  DEFAULT_MSG
  OPENCL_LIBRARY OPENCL_INCLUDE_DIR)

if(OPENCL_FOUND)
  set(OPENCL_LIBRARIES ${OPENCL_LIBRARY})
else(OPENCL_FOUND)
  set(OPENCL_LIBRARIES)
endif(OPENCL_FOUND)

mark_as_advanced(
  OPENCL_INCLUDE_DIR
  OPENCL_LIBRARY
  )
  
INCLUDE_DIRECTORIES(${OPENCL_INCLUDE_DIR})

include(CheckCompilerFlag)

CHECK_COMPILER_FLAG(-mavx AVX_FLAG)
CHECK_COMPILER_FLAG(-mavx2 AVX2_FLAG)

if (AVX_FLAG AND AVX2_FLAG)
	message("AVX2 support detected")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx -mavx2")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx -mavx2")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_CXX_FLAGS} -mavx -mavx2")
else()
	message("No AVX2 support detected")
endif()

add_subdirectory (util)
add_subdirectory (Kernels)
add_subdirectory (Timer)
add_subdirectory (impl)