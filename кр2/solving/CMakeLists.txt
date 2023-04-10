cmake_minimum_required(VERSION 3.25)

project(SLAE2_0)

set(CMAKE_CXX_STANDARD 20)

file(GLOB_RECURSE  solving_files *.h)
add_library(SLAE2.0 INTERFACE ${solving_files})