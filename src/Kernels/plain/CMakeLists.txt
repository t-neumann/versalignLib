MESSAGE( STATUS "CMAKE_CURRENT_BINARY_DIR: " ${CMAKE_CURRENT_BINARY_DIR} )
MESSAGE( STATUS "CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR} )

add_library(plainKernel SHARED SWKernel.cpp SWKernel_dllexport.cpp)
target_link_libraries (plainKernel Utils)
