MESSAGE( STATUS "CMAKE_CURRENT_BINARY_DIR: " ${CMAKE_CURRENT_BINARY_DIR} )
MESSAGE( STATUS "CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR} )

include_directories (${PROJECT_BINARY_DIR})

MESSAGE (STATUS "project binary dir" ${PROJECT_BINARY_DIR})
add_executable(hexify hexify/hexify.cpp)

ADD_CUSTOM_COMMAND(
	OUTPUT ${PROJECT_BINARY_DIR}/test_kernel.h
	COMMAND $<TARGET_FILE:hexify> ${CMAKE_CURRENT_SOURCE_DIR}/opencl_c/test_kernel
	COMMAND ${CMAKE_COMMAND} -E rename ${CMAKE_CURRENT_SOURCE_DIR}/opencl_c/test_kernel.h ${PROJECT_BINARY_DIR}/test_kernel.h
	DEPENDS hexify
	COMMENT "Hexifying test_kernel"
	MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/opencl_c/test_kernel.cl
)

add_library(openCLKernel STATIC ${PROJECT_BINARY_DIR}/test_kernel.h ocl_testing.cpp)	



