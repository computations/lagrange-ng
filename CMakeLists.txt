cmake_minimum_required(VERSION 3.13)
project(lagrange-cpp)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

set(LAGRANGE_SOURCES
        src/AncSplit.cpp
        src/Context.cpp
        src/InputReader.cpp
        src/Node.cpp
        src/Tree.cpp
        src/TreeReader.cpp
        src/Utils.cpp
        src/Workspace.cpp
        src/Operation.cpp
        src/WorkerState.cpp
        )

add_library(lagrange-release-libs OBJECT ${LAGRANGE_SOURCES})
add_executable(lagrange src/main.cpp)


find_package(PkgConfig)
pkg_search_module(MKL mkl-dynamic-lp64-iomp)

option(ENABLE_MKL "Enable MKL for linalg computations" ON)

if(${MKL_FOUND} AND ${ENABLE_MKL})
  set(LAGRANGE_LINALG_LDFLAGS ${MKL_LDFLAGS})
  set(LAGRANGE_LINALG_INCLUDE_DIRS ${MKL_INCLUDE_DIRS})
  set(LAGRANGE_LINALG_COMPILE_DEFS MKL_ENABLED) 
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

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -mtune=native -march=native -fno-omit-frame-pointer")
set(CMAKE_C_FLAGS_RELEASE "-O3 -mtune=native -march=native -fno-omit-frame-pointer")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -fno-omit-frame-pointer -fsanitize=address -fsanitize-recover=address")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g -fno-omit-frame-pointer -fsanitize=address -fsanitize-recover=address")
set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address -fsanitize-recover=address")

set(OPT_LIBS nlopt)
set(LAGRANGE_LIBS ${OPT_LIBS} nlohmann_json::nlohmann_json)

add_subdirectory(deps)

target_compile_options(lagrange-release-libs PRIVATE ${LAGRANGE_RELEASE_FLAGS} ${WARNING_FLAGS})
target_include_directories(lagrange-release-libs PRIVATE ${LAGRANGE_LINALG_INCLUDE_DIRS})
target_link_libraries(lagrange-release-libs PRIVATE ${LAGRANGE_LIBS} pthread
  ${LAGRANGE_LINALG_LDFLAGS})
target_compile_definitions(lagrange-release-libs PRIVATE
  ${LAGRANGE_LINALG_COMPILE_DEFS})

target_compile_options(lagrange PRIVATE ${LAGRANGE_RELEASE_FLAGS} ${WARNING_FLAGS})
target_include_directories(lagrange PRIVATE ${LAGRANGE_LINALG_INCLUDE_DIRS})
target_link_libraries(lagrange PRIVATE lagrange-release-libs pthread
  nlohmann_json::nlohmann_json ${LAGRANGE_LINALG_LDFLAGS})
target_compile_definitions(lagrange PRIVATE ${LAGRANGE_LINALG_COMPILE_DEFS})

add_subdirectory(tests)
