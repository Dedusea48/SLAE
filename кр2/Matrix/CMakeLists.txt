cmake_minimum_required(VERSION 3.25)

project(SLAE2_1)

set(CMAKE_CXX_STANDARD 20)

file(GLOB_RECURSE matrix_files *.h)
add_library(SLAE2_0 INTERFACE ${matrix_files})