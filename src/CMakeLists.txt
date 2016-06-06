include_directories (.)
include_directories ("${PROJECT_SOURCE_DIR}/include")

add_subdirectory (util)
add_subdirectory (Kernels)
add_subdirectory (Timer)

add_executable (Versalign main.cpp)
target_link_libraries (Versalign Kernels)
target_link_libraries (Versalign Timer)
target_link_libraries (Versalign Utils)
