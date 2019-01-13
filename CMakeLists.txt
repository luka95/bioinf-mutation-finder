cmake_minimum_required(VERSION 3.6)
project (MyProject)

# added -fopenmp
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

set(SOURCE_FILES Hirschberg.cpp)

add_executable (MyProject ${SOURCE_FILES})