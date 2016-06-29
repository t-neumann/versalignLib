MESSAGE( STATUS "CMAKE_CURRENT_BINARY_DIR: " ${CMAKE_CURRENT_BINARY_DIR} )
MESSAGE( STATUS "CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR} )

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
    
add_library(SSEKernel SSEKernel.cpp)
#add_library(AVXKernel AVXKernel.cpp)
target_link_libraries (SSEKernel Utils)
#target_link_libraries (AVXKernel Utils)
