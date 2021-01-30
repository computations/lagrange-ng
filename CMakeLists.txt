cmake_minimum_required(VERSION 3.13)
project(lagrange-cpp)
set(CMAKE_CXX_STANDARD 14)

add_subdirectory(deps)
set(INCLUDE_DIRS /usr/local/include src ${ARMADILLO_INCLUDES})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

link_directories(/usr/local/lib)

enable_language(Fortran)

set(LAGRANGE_SOURCES
        src/AncSplit.cpp
        src/AncSplit.h
        src/Context.cpp
        src/Context.h
        src/InputReader.cpp
        src/InputReader.h
        src/node.cpp
        src/node.h
        src/superdouble.cpp
        src/superdouble.h
        src/tree.cpp
        src/tree.h
        src/tree_reader.cpp
        src/tree_reader.h
        src/Utils.cpp
        src/Utils.h
        src/Workspace.h
        src/Workspace.cpp
        src/Operation.h
        src/Operation.cpp
        )


add_library(lagrange-release-libs OBJECT ${LAGRANGE_SOURCES})
add_library(lagrange-debug-libs OBJECT ${LAGRANGE_SOURCES})
add_library(lagrange-profile-libs OBJECT ${LAGRANGE_SOURCES})
add_executable(lagrange src/main.cpp)
add_executable(lagrange-debug src/main.cpp)
add_executable(lagrange-profile src/main.cpp)

# These flags are required for the expmkit code
set(CMAKE_Fortran_FLAGS_DEBUG "-fallow-argument-mismatch")
set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch")

set(WARNING_FLAGS -Wall -Wextra -Wshadow)

set(LAGRANGE_DEBUG_FLAGS -O0 -g -fno-omit-frame-pointer)
list(APPEND ${LAGRANGE_DEBUG_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG})
set(LAGRANGE_PROFILE_FLAGS -O3 -fno-omit-frame-pointer)
list(APPEND ${LAGRANGE_PROFILE_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE})
set(LAGRANGE_RELEASE_FLAGS -O3 -DBLAZE_USE_VECTORIZATION=1)
list(APPEND ${LAGRANGE_RELEASE_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE})

set(LINALG_LIBS gsl cblas lapack)
set(OPT_LIBS nlopt)
set(LAGRANGE_LIBS ${LINALG_LIBS} ${OPT_LIBS})

target_compile_options(lagrange-release-libs PRIVATE ${LAGRANGE_RELEASE_FLAGS}
        ${WARNING_FLAGS})
target_include_directories(lagrange-release-libs PRIVATE ${BLAZE_INCLUDES})
target_link_libraries(lagrange-release-libs ${LAGRANGE_LIBS})

target_link_libraries(lagrange lagrange-release-libs)
target_compile_options(lagrange PRIVATE ${LAGRANGE_RELEASE_FLAGS}
        ${WARNING_FLAGS})
target_include_directories(lagrange PRIVATE ${BLAZE_INCLUDES})

target_compile_options(lagrange-debug-libs PRIVATE ${LAGRANGE_DEBUG_FLAGS}
        ${WARNING_FLAGS})
target_include_directories(lagrange-debug-libs PRIVATE ${BLAZE_INCLUDES})
target_link_libraries(lagrange-debug-libs ${LAGRANGE_LIBS})

target_link_libraries(lagrange-debug lagrange-debug-libs)
target_compile_options(lagrange-debug PRIVATE ${LAGRANGE_DEBUG_FLAGS}
        ${WARNING_FLAGS})
target_include_directories(lagrange-debug PRIVATE ${BLAZE_INCLUDES})

target_compile_options(lagrange-profile-libs PRIVATE ${LAGRANGE_PROFILE_FLAGS}
        ${WARNING_FLAGS})
target_include_directories(lagrange-profile-libs PRIVATE ${BLAZE_INCLUDES})
target_link_libraries(lagrange-profile-libs ${LAGRANGE_LIBS})

target_link_libraries(lagrange-profile lagrange-profile-libs)
target_compile_options(lagrange-profile PRIVATE ${LAGRANGE_PROFILE_FLAGS}
        ${WARNING_FLAGS})
target_include_directories(lagrange-profile PRIVATE ${BLAZE_INCLUDES})

add_subdirectory(tests)
