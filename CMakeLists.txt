cmake_minimum_required(VERSION 3.24)
cmake_policy(SET CMP0135 NEW)

project(lagrange-ng)
set(CMAKE_CXX_STANDARD 23 REQUIRED)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")


include(FetchContent)

option(BUILD_SHARED_LIBS "Build using shared libraries" OFF)
option(BUILD_ONLY_STATIC "Only build the exclusively static library" OFF)

if("${CMAKE_BUILD_TYPE}" MATCHES "Debug" AND ${BUILD_ONLY_STATIC})
  message(WARNING "Cannot build a debug static binary")
endif()

# Find NLopt, otherwise download it
FetchContent_Declare(NLopt
  URL https://github.com/stevengj/nlopt/archive/refs/tags/v2.10.0.tar.gz
  URL_HASH SHA256=506f83a9e778ad4f204446e99509cb2bdf5539de8beccc260a014bd560237be1
  EXCLUDE_FROM_ALL
  FIND_PACKAGE_ARGS
)


# Find json lib, otherwise download it
FetchContent_Declare(nlohmann_json
  URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz
  URL_HASH SHA256=d6c65aca6b1ed68e7a182f4757257b107ae403032760ed6ef121c9d55e81757d
  EXCLUDE_FROM_ALL
  FIND_PACKAGE_ARGS
)

FetchContent_Declare(logger
  URL https://github.com/computations/logger/archive/refs/tags/v0.4.2.4.tar.gz
  URL_HASH SHA256=4ab8a9a044c682aa60107453567109137eff0c5e2c4c930ef715d2900baf3545
  EXCLUDE_FROM_ALL
  FIND_PACKAGE_ARGS
)

FetchContent_MakeAvailable(NLopt nlohmann_json logger)

FetchContent_GetProperties(NLopt
  POPULATED NLopt_POPULATED
)

if(NLopt_POPULATED)
  set(NLOPT_LIB nlopt)
else()
  set(NLOPT_LIB NLopt::nlopt)
endif()


# Find MKL, fallback to OpenBLAS
#include(FindMKL)
find_package(MKL CONFIG)

if (MKL_FOUND)
  set(LAGRANGE_LIBS $<LINK_ONLY:MKL::MKL>)
  set(LAGRANGE_INCLUDES $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
  set(LAGRANGE_BUILD_FLAGS $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
  set(LAGRANGE_DEFINES MKL_ENABLED)
else()
    if(DEFINED OpenBLAS_INCLUDE_DIRS)
      set(LAGRANGE_LIBS openblas)

      include(CheckSourceCompiles)
      include(CMakePushCheckState)

      cmake_push_check_state()

      set(CMAKE_REQUIRED_LIBRARIES openblas)
      set(CMAKE_REQUIRED_FLAGS -fopenmp)
      set(CMAKE_REQUIRED_INCLUDES ${OpenBLAS_INCLUDE_DIRS})

      check_source_compiles(CXX 
      "
      #include <cblas.h>
      #include <lapacke.h>

      int main(){
          LAPACKE_zgesv(CblasRowMajor,
                  1,
                  1,
                  nullptr,
                  1,
                  nullptr,
                  nullptr,
                  1);
          return 0;

      }
      "
      HAVE_LAPACKE)

      if(NOT HAVE_LAPACKE)

          set(CMAKE_REQUIRED_LIBRARIES openblas lapacke)

          check_source_compiles(CXX 
          "
          #include <cblas.h>
          #include <lapacke.h>

          int main(){
              LAPACKE_zgesv(CblasRowMajor,
                      1,
                      1,
                      nullptr,
                      1,
                      nullptr,
                      nullptr,
                      1);
              return 0;

          }
          "
          HAVE_LAPACKE_IN_SYSTEM_PATH)

        set(LAGRANGE_LIBS ${LAGRANGE_LIBS} lapacke)
      endif()

      cmake_pop_check_state()

    else()
      FetchContent_Declare(OpenBLAS
        URL
        https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.30/OpenBLAS-0.3.30.tar.gz
        URL_HASH SHA256=27342cff518646afb4c2b976d809102e368957974c250a25ccc965e53063c95d
        EXCLUDE_FROM_ALL
        FIND_PACKAGE_ARGS
      )

      FetchContent_MakeAvailable(OpenBLAS)
      set(OpenBLAS_INCLUDE_DIRS
        ${CMAKE_CURRENT_BINARY_DIR}/
        ${CMAKE_CURRENT_BINARY_DIR}/generated
        ${openblas_SOURCE_DIR}/lapack-netlib/LAPACKE/include/
      )

      message("Setting OpenBLAS library to ${OpenBLAS_LIBRARY_DIRS}")
      message("Setting OpenBLAS include directory to ${OpenBLAS_INCLUDE_DIRS}")
      set(LAGRANGE_LIBS openblas_static)
      endif()

      set(LAGRANGE_LIB_DIR ${OpenBLAS_LIBRARY_DIRS})
      set(LAGRANGE_INCLUDES ${OpenBLAS_INCLUDE_DIRS})
endif()

find_package(Git)
if(GIT_FOUND)
    if(EXISTS ${CMAKE_SOURCE_DIR})
        execute_process(
            COMMAND ${GIT_EXECUTABLE} describe
            WORKING_DIRECTORY "${local_dir}"
            OUTPUT_VARIABLE LAGRANGE_BUILD_VERSION
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE)
    endif()
else()
  SET(LAGRANGE_BUILD_VERSION "v0.7.0-release")
endif()

set(LAGRANGE_DEFINES ${LAGRANGE_DEFINES} LAGRANGE_BUILD_VERSION=${LAGRANGE_BUILD_VERSION})

add_subdirectory(src)

option(BUILD_TESTS "Enable building the tests" OFF)

if(BUILD_TESTS AND NOT BUILD_ONLY_STATIC)
  include(CTest)
  enable_testing()
  add_subdirectory(tests)
endif()
