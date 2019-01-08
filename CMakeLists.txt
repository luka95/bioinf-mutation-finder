cmake_minimum_required(VERSION 3.13)
project(bioinf_mutation_finder)

set(CMAKE_CXX_STANDARD 17)

add_executable(bioinf_mutation_finder   mapping_demo.cpp
                                        Index.cpp Index.h
                                        LIS.cpp LIS.h
                                        Align.cpp Align.h
                                        DataFetcher.cpp
                                        DataFetcher.h
                                        Mutation.cpp
                                        Mutation.h
                                        tests/index_test.cpp
                                        tests/LIS_test.cpp)
add_subdirectory(tests)