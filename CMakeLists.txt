cmake_minimum_required(VERSION 3.12)
project(Bioinformatika)

set(CMAKE_CXX_STANDARD 14)

add_executable(Bioinformatika main.cpp Mutation.cpp Mutation.h DataFetcher.cpp DataFetcher.h Base.h Functions.cpp Functions.h Align.h Align.cpp Index.cpp Index.h LIS.h LIS.cpp)