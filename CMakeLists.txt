cmake_minimum_required(VERSION 3.5)

project(TP_01_Jonathan_Guimaraes_openmp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(OpenMP)

if (OPENMP_FOUND)
    message("Openmp Founded")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_subdirectory(nbody_parallel)
add_executable(TP_01_Jonathan_Guimaraes_openmp test.cpp)

target_link_libraries(TP_01_Jonathan_Guimaraes_openmp PUBLIC nbody_parallel)
target_include_directories(TP_01_Jonathan_Guimaraes_openmp PUBLIC "${CMAKE_BUILD_RPATH}/nbody_parallel/")
