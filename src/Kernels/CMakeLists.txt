add_library(Kernels SWKernel.cpp SSEKernel.cpp)
target_link_libraries (Kernels Utils)
