cmake_minimum_required(VERSION 2.10)

project(mqlib)

set(CMAKE_CXX_STANDARD 11)

set(SOURCES
        src/util/randomForest.cpp
        src/util/random.cpp
        src/problem/qubo_instance.cpp
        src/problem/qubo_heuristic.cpp
        src/problem/max_cut_instance.cpp
        src/problem/max_cut_heuristic.cpp
        src/problem/instance.cpp
        src/problem/heuristic.cpp
        src/metrics/max_cut_metrics.cpp
        src/heuristics/heuristic_factory.cpp
        src/heuristics/extended_solution.cpp
        src/heuristics/base_solution.cpp
        src/heuristics/qubo/qubo_solution.cpp
        src/heuristics/qubo/qubo_simple_solution.cpp
        src/heuristics/qubo/qubo_partial_solution.cpp
        src/heuristics/qubo/pardalos2008.cpp
        src/heuristics/qubo/palubeckis2006.cpp
        src/heuristics/qubo/palubeckis2004b.cpp
        src/heuristics/qubo/merz2004.cpp
        src/heuristics/qubo/merz2002.cpp
        src/heuristics/qubo/merz1999.cpp
        src/heuristics/qubo/lu2010.cpp
        src/heuristics/qubo/lodi1999.cpp
        src/heuristics/qubo/katayama2001.cpp
        src/heuristics/qubo/katayama2000.cpp
        src/heuristics/qubo/hasan2000.cpp
        src/heuristics/qubo/glover2010.cpp
        src/heuristics/qubo/glover1998a.cpp
        src/heuristics/qubo/beasley1998.cpp
        src/heuristics/qubo/alkhamis1998.cpp
        src/heuristics/maxcut/max_cut_solution.cpp
        src/heuristics/maxcut/max_cut_simple_solution.cpp
        src/heuristics/maxcut/max_cut_partial_solution.cpp
        src/heuristics/maxcut/laguna2009.cpp
        src/heuristics/maxcut/hyperheuristic.cpp
        src/heuristics/maxcut/festa2002.cpp
        src/heuristics/maxcut/duarte2005.cpp
        src/heuristics/maxcut/deSousa2013.cpp
        src/heuristics/maxcut/burer2002.cpp
        src/heuristics/maxcut/baseline.cpp
        )

add_library(mqlib SHARED ${SOURCES})

target_include_directories(mqlib PUBLIC ${PROJECT_SOURCE_DIR}/include)


# Build binary
option(MQLIB_BINARY "Build a cli binary" OFF)
if (MQLIB_BINARY)
    add_subdirectory(src)
endif ()

# Enable testing when requested by the user
option(MQLIB_TESTS "Build the tests" OFF)
if (MQLIB_TESTS)
    enable_testing()
    add_subdirectory(test)
endif ()