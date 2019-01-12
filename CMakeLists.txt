cmake_minimum_required(VERSION 3.13)
project(bioinf_mutation_finder)

set(CMAKE_CXX_STANDARD 14)

include_directories(src)
include_directories(src/alignment)
include_directories(src/mapping/index)
include_directories(src/mapping/inverter)
include_directories(src/mapping/lis)

add_executable(bioinf_mutation_finder
        src/alignment/Align.cpp
        src/alignment/Align.h
        src/alignment/align_demo.cpp
        src/alignment/Hirschberg.cpp
        src/mapping/index/Index.cpp
        src/mapping/index/Index.h
        src/mapping/index/index_demo.cpp
        src/mapping/inverter/Inverter.cpp
        src/mapping/inverter/Inverter.h
        src/mapping/inverter/inverter_demo.cpp
        src/mapping/lis/LIS.cpp
        src/mapping/lis/LIS.h
        src/mapping/lis/LIS_demo.cpp
        src/mapping/mapping_demo.cpp
        src/DataLoader.cpp
        src/DataLoader.h
        src/main.cpp
        src/Mutation.cpp
        src/Mutation.h)
