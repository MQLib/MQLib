file(GLOB_RECURSE source_files *.cpp)
list(REMOVE_ITEM source_files "main.cpp")

add_library(MQLib STATIC ${source_files})

add_executable( MQLib_bin main.cpp ${headers} )
target_link_libraries( MQLib_bin MQLib m )
set_target_properties(MQLib_bin PROPERTIES OUTPUT_NAME MQLib)
