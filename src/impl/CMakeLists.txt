cmake_minimum_required (VERSION 2.8.11)

add_executable (Versalign main.cpp)

target_link_libraries (Versalign Utils)
target_link_libraries (Versalign "-ldl")