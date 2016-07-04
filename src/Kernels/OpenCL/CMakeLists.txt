MESSAGE( STATUS "CMAKE_CURRENT_BINARY_DIR: " ${CMAKE_CURRENT_BINARY_DIR} )
MESSAGE( STATUS "CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR} )

#set(Boost_USE_STATIC_LIBS OFF) 
#set(Boost_USE_MULTITHREADED ON)  
#set(Boost_USE_STATIC_RUNTIME OFF) 
#find_package(Boost REQUIRED COMPONENTS compute) 

#if(Boost_FOUND)
#    include_directories(${Boost_INCLUDE_DIRS})
    add_library(openCLKernel ocl_testing.cpp) 
#    target_link_libraries(openCLKernel ${Boost_LIBRARIES})
#endif()

add_executable(hexify hexify/hexify.cpp)



