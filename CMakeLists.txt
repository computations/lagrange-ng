cmake_minimum_required(VERSION 3.13)
project(lagrange-cpp)
set(CMAKE_CXX_STANDARD 14)

add_subdirectory(deps)
set(INCLUDE_DIRS /usr/local/include src ${ARMADILLO_INCLUDES})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

link_directories(/usr/local/lib)

set(LAGRANGE_SOURCES
        src/AncSplit.cpp
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

set(LAGRANGE_RELEASE_FLAGS -DBLAZE_USE_VECTORIZATION=1)
list(APPEND ${LAGRANGE_RELEASE_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE})

#set(LINALG_LIBS gsl cblas lapack)
set(OPT_LIBS nlopt)
set(LAGRANGE_LIBS ${OPT_LIBS} nlohmann_json::nlohmann_json)


add_executable(lagrange src/main.cpp)
add_library(lagrange-release-libs OBJECT ${LAGRANGE_SOURCES})

add_dependencies(lagrange-release-libs OpenBLAS)
add_dependencies(lagrange OpenBLAS)

target_compile_options(lagrange-release-libs PRIVATE ${LAGRANGE_RELEASE_FLAGS}
    ${WARNING_FLAGS})
target_include_directories(lagrange-release-libs PRIVATE ${BLAZE_INCLUDES}
    ${OpenBLAS_INCLUDE_DIRS})
target_link_libraries(lagrange-release-libs PRIVATE ${LAGRANGE_LIBS}
    ${OpenBLAS_LIB} pthread)

target_link_libraries(lagrange PRIVATE lagrange-release-libs pthread
    nlohmann_json::nlohmann_json)
target_compile_options(lagrange PRIVATE ${LAGRANGE_RELEASE_FLAGS}
    ${WARNING_FLAGS})
target_include_directories(lagrange PRIVATE ${BLAZE_INCLUDES}
    ${OpenBLAS_INCLUDE_DIRS})


add_subdirectory(tests)
