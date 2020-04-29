cmake_minimum_required(VERSION 3.10)
project(SnowSim)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories(EIGEN3_INCLUDE_DIR)
include_directories(${CMAKE_SOURCE_DIR}/include/)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -pthread")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

add_executable(SnowSim
               src/main.cpp)