cmake_minimum_required(VERSION 2.6)
project(tests)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fopenmp")

set(SOURCE_FILES main.cpp TestMatrixMarketWriter.cpp ../output/MatrixMarketWriter.cpp)
add_executable(tests ${SOURCE_FILES})
