MESSAGE( STATUS "CMAKE_CURRENT_BINARY_DIR: " ${CMAKE_CURRENT_BINARY_DIR} )
MESSAGE( STATUS "CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR} )

include_directories (${PROJECT_BINARY_DIR})

MESSAGE (STATUS "project binary dir" ${PROJECT_BINARY_DIR})
add_executable(hexify hexify/hexify.cpp)

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

add_library(openCLKernel STATIC ${PROJECT_BINARY_DIR}/opencl_definitions.h ${PROJECT_BINARY_DIR}/scoring_kernels.h ${PROJECT_BINARY_DIR}/alignment_kernels.h ocl_testing.cpp)	



