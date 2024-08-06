cmake_minimum_required(VERSION 3.16)

project(dcpbwt)
add_executable(${PROJECT_NAME} main.cpp)

# ---- Download submodules ---- #
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()


add_subdirectory(external/DYNAMIC)
target_include_directories(${PROJECT_NAME} PUBLIC external/DYNAMIC/include)
target_link_libraries(${PROJECT_NAME} PRIVATE dynamic)
include_directories(${hopscotch_map_INCLUDE})
add_dependencies(${PROJECT_NAME} hopscotch_map)

#include_directories(${PROJECT_SOURCE_DIR}/build/external/DYNAMIC/hopscotch_map-prefix/src/hopscotch_map/include)
#message(PROJECT_SOURCE_DIR="${hopscotch_map_INCLUDE}")
