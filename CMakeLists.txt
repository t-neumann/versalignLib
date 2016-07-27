# CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8.11)
project (Versalign)

MESSAGE( STATUS "CMAKE_SYSTEM: " ${CMAKE_SYSTEM} )

set( VERSION_MAJOR 0 )
set( VERSION_MINOR 1 )

IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Building in debug mode.")
    set( VERSION_BUILD 0-debug )
else()
  set( VERSION_BUILD 0 )
ENDIF()


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin/Versalign-${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_BUILD}/)
file(MAKE_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})

#set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin/)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin/Versalign-${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_BUILD}/lib/)
file(MAKE_DIRECTORY ${LIBRARY_OUTPUT_PATH})
#set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin/)

# Recurse into the "Hello" and "Demo" subdirectories. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.
add_subdirectory (src)