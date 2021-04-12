cmake_minimum_required(VERSION 3.13)
project(lagrange-cpp)
set(CMAKE_CXX_STANDARD 14)

set(INCLUDE_DIRS /usr/local/include src ${ARMADILLO_INCLUDES})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

link_directories(/usr/local/lib)
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

set(LAGRANGE_RELEASE_FLAGS -DBLAZE_USE_VECTORIZATION)
list(APPEND ${LAGRANGE_RELEASE_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE})

#set(LINALG_LIBS gsl cblas lapack)
set(OPT_LIBS nlopt)
set(LAGRANGE_LIBS ${OPT_LIBS} nlohmann_json::nlohmann_json lapack)


add_executable(lagrange src/main.cpp)
add_library(lagrange-release-libs OBJECT ${LAGRANGE_SOURCES})

add_subdirectory(deps)

target_compile_options(lagrange-release-libs PRIVATE ${LAGRANGE_RELEASE_FLAGS}
    ${WARNING_FLAGS})
target_include_directories(lagrange-release-libs PRIVATE ${BLAZE_INCLUDES}
    ${BLAS_LIBS_INCLUDE_DIRS})
target_link_libraries(lagrange-release-libs PRIVATE ${LAGRANGE_LIBS}
    ${BLAS_LIBS_BINARY} pthread OpenMP::OpenMP_CXX lapacke)

target_link_libraries(lagrange PRIVATE lagrange-release-libs pthread
    nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX ${BLAS_LIBS_BINARY} lapacke)
target_compile_options(lagrange PRIVATE ${LAGRANGE_RELEASE_FLAGS} 
    ${WARNING_FLAGS})
target_include_directories(lagrange PRIVATE ${BLAZE_INCLUDES}
    ${BLAS_LIBS_INCLUDE_DIRS} )


add_subdirectory(tests)
