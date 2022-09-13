add_definitions(-DMQLIB_PROJECT_DIR="${CMAKE_SOURCE_DIR}")

add_executable(test_sample_problems test_sample_problems.cpp)
target_link_libraries(test_sample_problems mqlib)
