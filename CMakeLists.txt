cmake_minimum_required(VERSION 3.13)
project(lagrange-ng)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

set(LAGRANGE_SOURCES
        src/AncSplit.cpp
        src/Context.cpp
        src/Node.cpp
        src/Tree.cpp
        src/TreeReader.cpp
        src/Utils.cpp
        src/Workspace.cpp
        src/Operation.cpp
        src/WorkerState.cpp
        src/Arnoldi.cpp
        src/Fossil.cpp
        src/ConfigFile.cpp
        src/Alignment.cpp
        )

add_library(lagrange-release-libs OBJECT ${LAGRANGE_SOURCES})
add_executable(lagrange src/main.cpp)


find_package(PkgConfig)
pkg_search_module(MKL mkl-dynamic-lp64-iomp)

option(ENABLE_MKL "Enable MKL for linalg computations" ON)

set(MKL_PREFIX_DIR CACHE STRING "")

if((${ENABLE_MKL} AND ${MKL_FOUND}) OR (NOT "${MKL_PREFIX_DIR}" STREQUAL ""))
  if(${MKL_FOUND})
    set(LAGRANGE_LINALG_LDFLAGS ${MKL_LDFLAGS})
    set(LAGRANGE_LINALG_INCLUDE_DIRS ${MKL_INCLUDE_DIRS})
    set(LAGRANGE_LINALG_COMPILE_DEFS MKL_ENABLED) 
  else()
    set(LAGRANGE_LINALG_LDFLAGS mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 m dl)
    set(LAGRANGE_LINALG_INCLUDE_DIRS ${MKL_PREFIX_DIR}/include/)
    set(LAGRANGE_LINALG_LINK_DIRS ${MKL_PREFIX_DIR}/lib/)
    set(LAGRANGE_LINALG_COMPILE_DEFS MKL_ENABLED) 
  endif()
else()
  find_package(OpenMP REQUIRED)
  include(ExternalProject)

  set(OpenBLAS_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/OpenBLAS)

  ExternalProject_Add(OpenBLAS
  GIT_REPOSITORY https://github.com/xianyi/OpenBLAS
  GIT_TAG v0.3.19
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE
  CONFIGURE_COMMAND ""
  BUILD_IN_SOURCE TRUE
  BUILD_COMMAND make USE_OPENMP=1
  INSTALL_COMMAND make PREFIX=${OpenBLAS_INSTALL_PREFIX} install
  UPDATE_DISCONNECTED TRUE
  )
  add_dependencies(lagrange-release-libs OpenBLAS)
  add_dependencies(lagrange OpenBLAS)
  
  set(LAGRANGE_LINALG_LDFLAGS ${OpenBLAS_INSTALL_PREFIX}/lib/libopenblas.a
    OpenMP::OpenMP_CXX)
  set(LAGRANGE_LINALG_INCLUDE_DIRS ${OpenBLAS_INSTALL_PREFIX}/include)
  set(LAGRANGE_LINALG_COMPILE_DEFS) 
endif()

set(WARNING_FLAGS -Wall -Wextra -Wshadow)

if(${CMAKE_BUILD_TYPE} AND ${CMAKE_BUILD_TYPE} EQUAL "Release")
  set(LAGRANGE_BUILD_FLAGS "-O3 -mtune=native -march=native -fno-omit-frame-pointer")
elseif(${CMAKE_BUILD_TYPE} AND ${CMAKE_BUILD_TYPE} EQUAL "Debug")
  set(LAGRANGE_BUILD_FLAGS "-O0 -g -fno-omit-frame-pointer -fsanitize=address -fsanitize-recover=address")
endif()

set(OPT_LIBS nlopt)
set(LAGRANGE_LIBS ${OPT_LIBS} nlohmann_json::nlohmann_json)

add_subdirectory(deps)

target_compile_options(lagrange-release-libs PUBLIC ${LAGRANGE_BUILD_FLAGS} ${WARNING_FLAGS})
target_include_directories(lagrange-release-libs PUBLIC ${LAGRANGE_LINALG_INCLUDE_DIRS})
target_link_libraries(lagrange-release-libs PUBLIC ${LAGRANGE_LIBS} pthread
  ${LAGRANGE_LINALG_LDFLAGS})
target_compile_definitions(lagrange-release-libs PUBLIC
  ${LAGRANGE_LINALG_COMPILE_DEFS})
target_link_directories(lagrange-release-libs PUBLIC ${LAGRANGE_LINALG_LINK_DIRS})

target_compile_options(lagrange PRIVATE ${LAGRANGE_BUILD_FLAGS} ${WARNING_FLAGS})
target_include_directories(lagrange PRIVATE ${LAGRANGE_LINALG_INCLUDE_DIRS})
target_link_libraries(lagrange PRIVATE lagrange-release-libs pthread
  nlohmann_json::nlohmann_json ${LAGRANGE_LINALG_LDFLAGS})
target_compile_definitions(lagrange PRIVATE ${LAGRANGE_LINALG_COMPILE_DEFS})
target_link_directories(lagrange PRIVATE ${LAGRANGE_LINALG_LINK_DIRS})

add_subdirectory(tests)
