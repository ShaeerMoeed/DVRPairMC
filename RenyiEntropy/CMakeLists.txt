cmake_minimum_required(VERSION 3.27)
project(RenyiEntropy)

set(CMAKE_CXX_STANDARD 20)

add_executable(RenyiEntropy main.cpp
        OutputWriter.cpp Grid.cpp)
