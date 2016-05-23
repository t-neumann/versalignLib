include_directories (.)
include_directories ("${PROJECT_SOURCE_DIR}/include")

add_subdirectory (Kernels)

add_executable (Versalign main.cpp)
target_link_libraries (Versalign SWKernel)