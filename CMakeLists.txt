cmake_minimum_required(VERSION 3.10)
project(SnowSim)

#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

include_directories(EIGEN3_INCLUDE_DIR)
include_directories(${CMAKE_SOURCE_DIR}/include/)
include_directories(${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS})

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -pthread")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

add_executable(SnowSim
               src/main.cpp
               src/BVH.cpp
               src/camera.cpp
               src/SnowParticle.cpp)

target_link_libraries(SnowSim ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})