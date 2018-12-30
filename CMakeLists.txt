cmake_minimum_required(VERSION 3.13)
project(bioinf_mutation_finder)

set(CMAKE_CXX_STANDARD 17)

add_executable(bioinf_mutation_finder index_demo.cpp Index.cpp Index.h)