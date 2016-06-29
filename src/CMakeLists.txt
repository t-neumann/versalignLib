MESSAGE( STATUS "CMAKE_CURRENT_BINARY_DIR: " ${CMAKE_CURRENT_BINARY_DIR} )
MESSAGE( STATUS "CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR} )

include_directories (.)
include_directories ("${PROJECT_SOURCE_DIR}/include")

add_subdirectory (util)
add_subdirectory (Kernels)
add_subdirectory (Timer)

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

add_executable (Versalign main.cpp)
target_link_libraries (Versalign plainKernel)
target_link_libraries (Versalign openCLKernel)
target_link_libraries (Versalign SSEKernel)
target_link_libraries (Versalign Timer)
target_link_libraries (Versalign Utils)

target_link_libraries(Versalign ${OPENCL_LIBRARY})
