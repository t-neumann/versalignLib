cmake_minimum_required (VERSION 2.8.11)

include_directories (${PROJECT_BINARY_DIR})

###################################################
# Build hexify executable for OpenCL source import
###################################################

MESSAGE (STATUS "project binary dir" ${PROJECT_BINARY_DIR})
add_executable(hexify hexify/hexify.cpp)

###################################################
# Import OpenCL source
###################################################

ADD_CUSTOM_COMMAND(
	OUTPUT ${PROJECT_BINARY_DIR}/opencl_definitions.h
	COMMAND $<TARGET_FILE:hexify> ${CMAKE_CURRENT_SOURCE_DIR}/opencl_c/opencl_definitions
	COMMAND ${CMAKE_COMMAND} -E rename ${CMAKE_CURRENT_SOURCE_DIR}/opencl_c/opencl_definitions.h ${PROJECT_BINARY_DIR}/opencl_definitions.h
	DEPENDS hexify
	COMMENT "Hexifying opencl_definitions"
	MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/opencl_c/opencl_definitions.cl
)

ADD_CUSTOM_COMMAND(
	OUTPUT ${PROJECT_BINARY_DIR}/scoring_kernels.h
	COMMAND $<TARGET_FILE:hexify> ${CMAKE_CURRENT_SOURCE_DIR}/opencl_c/scoring_kernels
	COMMAND ${CMAKE_COMMAND} -E rename ${CMAKE_CURRENT_SOURCE_DIR}/opencl_c/scoring_kernels.h ${PROJECT_BINARY_DIR}/scoring_kernels.h
	DEPENDS hexify
	COMMENT "Hexifying scoring_kernels"
	MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/opencl_c/scoring_kernels.cl
)

ADD_CUSTOM_COMMAND(
	OUTPUT ${PROJECT_BINARY_DIR}/alignment_kernels.h
	COMMAND $<TARGET_FILE:hexify> ${CMAKE_CURRENT_SOURCE_DIR}/opencl_c/alignment_kernels
	COMMAND ${CMAKE_COMMAND} -E rename ${CMAKE_CURRENT_SOURCE_DIR}/opencl_c/alignment_kernels.h ${PROJECT_BINARY_DIR}/alignment_kernels.h
	DEPENDS hexify
	COMMENT "Hexifying alignment_kernels"
	MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/opencl_c/alignment_kernels.cl
)


###################################################
# Link OpenCL libraries
###################################################

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

add_library(OpenCLKernel SHARED ${PROJECT_BINARY_DIR}/opencl_definitions.h ${PROJECT_BINARY_DIR}/scoring_kernels.h ${PROJECT_BINARY_DIR}/alignment_kernels.h OpenCLKernel.cpp OpenCLKernel_dllexport.cpp)
target_link_libraries(OpenCLKernel ${OPENCL_LIBRARY})