cmake_minimum_required(VERSION 3.16)

project(dcpbwt_new)
add_executable(${PROJECT_NAME} main.cpp)


add_subdirectory(external/DYNAMIC)

target_include_directories(${PROJECT_NAME} PUBLIC external/DYNAMIC/include)
target_link_libraries(${PROJECT_NAME} PRIVATE dynamic)

#include_directories(${PROJECT_SOURCE_DIR}/build/external/DYNAMIC/hopscotch_map-prefix/src/hopscotch_map/include)
message(PROJECT_SOURCE_DIR="${hopscotch_map_INCLUDE}")

include_directories(${hopscotch_map_INCLUDE})

#target_link_directories()

add_dependencies(${PROJECT_NAME} hopscotch_map)
