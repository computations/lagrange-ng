cmake_minimum_required(VERSION 3.13)
project(lagrange-ng)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

include(FindMKL)

find_package(NLopt REQUIRED)
find_package(MKL REQUIRED)

add_subdirectory(lib)
add_subdirectory(src)

option(BUILD_TESTS "Enable building the tests" OFF)

if(${BUILD_TESTS})
  add_subdirectory(tests)
endif()
