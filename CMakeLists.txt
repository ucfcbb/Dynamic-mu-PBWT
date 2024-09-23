cmake_minimum_required(VERSION 3.16)
project(dcpbwt)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(OpenMP)
set(EXE_NAME dcpbwt)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  message("Using ${CMAKE_BUILD_TYPE} mode ")
  add_compile_options(-Wall -Wextra -pedantic -O0)
else()
  message("Using ${CMAKE_BUILD_TYPE} mode ")
  add_compile_options(-Wall -Wextra -O3 -march=native -fopenmp -ffast-math ${OpenMP_CXX_FLAGS})
endif()

add_executable(${EXE_NAME} main.cpp include/dcpbwt.h
        include/phi.h
        include/dcpbwt_column.h
)
add_executable(test_insert test_insertion.cpp include/dcpbwt.h
        include/phi.h
        include/dcpbwt_column.h)
add_executable(test_delete test_deletion.cpp include/dcpbwt.h
        include/phi.h
        include/dcpbwt_column.h)

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

# External Libraries
# DYNAMIC
add_subdirectory(external/DYNAMIC)

# Build htslib
include(FetchContent)
FetchContent_Declare(
        htslib
        GIT_REPOSITORY https://github.com/samtools/htslib.git
)
FetchContent_GetProperties(htslib)
if(NOT htslib_POPULATED)
    set(disable_flags --disable-gcs --disable-s3 --disable-plugins --disable-lzma
            --disable-bz2 --disable-libcurl --without-libdeflate)
    FetchContent_Populate(htslib)
    message("HTSLIB_SRC: ${htslib_SOURCE_DIR}")
    message("HTSLIB_BIN: ${htslib_BINARY_DIR}")
    execute_process(COMMAND autoreconf -i
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${htslib_SOURCE_DIR})
    if (result)
        message(FATAL_ERROR "autoreconf step for htslib failed: ${result}")
    endif ()
    execute_process(COMMAND ./configure --prefix=${CMAKE_CURRENT_BINARY_DIR} ${disable_flags}
            WORKING_DIRECTORY ${htslib_SOURCE_DIR})
    if (result)
        message(FATAL_ERROR "configure step for htslib failed: ${result}")
    endif ()
    execute_process(COMMAND make lib-static
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${htslib_SOURCE_DIR})
    if (result)
        message(FATAL_ERROR "make step for htslib failed: ${result}")
    endif ()
    execute_process(COMMAND make install prefix=${CMAKE_CURRENT_BINARY_DIR}
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${htslib_SOURCE_DIR})
    if (result)
        message(FATAL_ERROR "install step for htslib failed: ${result}")
    endif()
    set(htslib_SRC ${CMAKE_CURRENT_BINARY_DIR}/include)
    set(htslib_LIB ${CMAKE_CURRENT_BINARY_DIR}/lib/libhts.a)
endif()
add_library(htslib STATIC IMPORTED GLOBAL)
set_property(TARGET htslib PROPERTY IMPORTED_LOCATION ${htslib_LIB})
target_include_directories(htslib INTERFACE ${htslib_SRC})

target_include_directories(${EXE_NAME} PUBLIC external/DYNAMIC/include PUBLIC include)
target_include_directories(test_insert PUBLIC external/DYNAMIC/include PUBLIC include)
target_include_directories(test_delete PUBLIC external/DYNAMIC/include PUBLIC include)

# DYNAMIC's hopscotch dependency
include_directories(${CMAKE_CURRENT_BINARY_DIR}/external/DYNAMIC/hopscotch_map-prefix/src/hopscotch_map/include)

target_link_libraries(${EXE_NAME} ${OpenMP_CXX_LIBRARIES})
target_link_libraries(${EXE_NAME} htslib -lz)
target_link_libraries(test_insert htslib -lz)
target_link_libraries(test_delete htslib -lz)

add_dependencies(${EXE_NAME} hopscotch_map)
add_dependencies(test_insert hopscotch_map)
add_dependencies(test_delete hopscotch_map)

