cmake_minimum_required(VERSION 3.25)
project(SLAE2_0)

set(CMAKE_CXX_STANDARD 20)

enable_testing()

add_executable(SLAE KR2.cpp)
add_subdirectory(Matrix)
add_subdirectory(solving)
target_link_libraries(SLAE2_0)