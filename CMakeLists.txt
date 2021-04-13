cmake_minimum_required(VERSION 3.13)
project(lagrange-cpp)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(OpenMP)

set(LAGRANGE_SOURCES
        src/AncSplit.cpp
        src/BranchSegment.cpp
        src/Context.cpp
        src/InputReader.cpp
        src/node.cpp
        src/tree.cpp
        src/tree_reader.cpp
        src/Utils.cpp
        src/Workspace.cpp
        src/Operation.cpp
        )

set(WARNING_FLAGS -Wall -Wextra -Wshadow)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fno-omit-frame-pointer")
set(CMAKE_C_FLAGS_RELEASE "-O3 -fno-omit-frame-pointer")

set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")

add_subdirectory(deps)

add_library(lagrange-release-libs OBJECT ${LAGRANGE_SOURCES})
add_executable(lagrange src/main.cpp)

target_compile_options(lagrange-release-libs PRIVATE ${WARNING_FLAGS})
target_include_directories(lagrange-release-libs PRIVATE)
target_link_libraries(lagrange-release-libs PRIVATE nlopt Eigen3::Eigen 
    OpenMP::OpenMP_CXX nlohmann_json::nlohmann_json)

target_compile_options(lagrange PRIVATE ${WARNING_FLAGS})
target_include_directories(lagrange PRIVATE)
target_link_libraries(lagrange PRIVATE lagrange-release-libs 
    nlopt Eigen3::Eigen OpenMP::OpenMP_CXX nlohmann_json::nlohmann_json)


add_subdirectory(tests)
