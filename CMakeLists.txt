project(MQLib)
cmake_minimum_required(VERSION 2.10)

set(MQLib_VERSION_MAJOR 0)
set(MQLib_VERSION_MINOR 1)

# C++11
add_compile_options(-std=c++11)

include_directories(include)
file(GLOB_RECURSE headers include/*.h)
add_subdirectory(src)
add_subdirectory(test)
