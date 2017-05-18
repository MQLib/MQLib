cmake_minimum_required(VERSION 2.8)
project(MQLib)

set(MQLib_VERSION_MAJOR 0)
set(MQLib_VERSION_MINOR 1)

# C++11
add_compile_options(-std=c++11)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")

include_directories(include)
file(GLOB_RECURSE headers include/*.h)
add_subdirectory(src)
