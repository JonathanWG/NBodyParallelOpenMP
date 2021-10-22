cmake_minimum_required(VERSION 3.5)

project(TP_01_Jonathan_Guimaraes_openmp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_subdirectory(nbody_parallel)
add_executable(TP_01_Jonathan_Guimaraes_openmp test.cpp)

target_link_libraries(TP_01_Jonathan_Guimaraes_openmp PUBLIC nbody_parallel)
target_include_directories(TP_01_Jonathan_Guimaraes_openmp PUBLIC "${CMAKE_BUILD_RPATH}/nbody_parallel/")
