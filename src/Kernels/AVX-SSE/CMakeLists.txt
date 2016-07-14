MESSAGE( STATUS "CMAKE_CURRENT_BINARY_DIR: " ${CMAKE_CURRENT_BINARY_DIR} )
MESSAGE( STATUS "CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR} )

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx -mavx2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx -mavx2")
    
add_library(SSEKernel SHARED SSEKernel.cpp SSEKernel_dllexport.cpp)
add_library(AVXKernel SHARED AVXKernel.cpp AVXKernel_dllexport.cpp)