file(GLOB_RECURSE source_files *.cpp)
list(REMOVE_ITEM source_files "main.cpp")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")

add_library(MQLib STATIC ${source_files})

add_executable(MQLib_bin main.cpp ${headers})
target_link_libraries(MQLib_bin MQLib m)
set_target_properties(MQLib_bin PROPERTIES OUTPUT_NAME MQLib)
