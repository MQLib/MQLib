set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

add_executable(MQLib_bin main.cpp)
target_link_libraries(MQLib_bin mqlib)
set_target_properties(MQLib_bin PROPERTIES OUTPUT_NAME MQLib)
