cmake_minimum_required(VERSION 3.24)
cmake_policy(SET CMP0135 NEW)

project(lagrange-ng)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")


include(FetchContent)

option(BUILD_SHARED_LIBS "Build using shared libraries" OFF)

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
  URL https://github.com/computations/logger/archive/refs/tags/v0.4.2.2.tar.gz
  URL_HASH SHA256=771e56674b73d013825273c300e19417f3b33944ce94589529cca3ebce581c7f
  EXCLUDE_FROM_ALL
  FIND_PACKAGE_ARGS
)

FetchContent_MakeAvailable(NLopt nlohmann_json logger)

FetchContent_GetProperties(NLopt
  POPULATED NLopt_POPULATED
)

if(${NLopt_POPULATED})
  set(NLOPT_LIB nlopt)
else()
  set(NLOPT_LIB NLopt::nlopt)
endif()


# Find MKL, fallback to OpenBLAS
include(FindMKL)
find_package(MKL)

if (${MKL_FOUND})
  set(LAGRANGE_LIBS ${MKL_LIBRARIES})
  set(LAGRANGE_INCLUDES ${MKL_INCLUDE_DIRS})
  set(LAGRANGE_DEFINES MKL_ENABLED=1)
else()
    FetchContent_Declare(OpenBLAS
      URL https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.29/OpenBLAS-0.3.29.tar.gz
      URL_HASH SHA256=38240eee1b29e2bde47ebb5d61160207dc68668a54cac62c076bb5032013b1eb
      EXCLUDE_FROM_ALL
      FIND_PACKAGE_ARGS
    )
    FetchContent_MakeAvailable(OpenBLAS)
    set(LAGRANGE_LIBS openblas)

    FetchContent_GetProperties(OpenBLAS
      POPULATED OpenBLAS_POPULATED
    )

    if(NOT ${OpenBLAS_FOUND})
      set(OpenBLAS_INCLUDE_DIRS 
        ${CMAKE_BINARY_DIR}/
        ${CMAKE_CURRENT_BINARY_DIR}/generated/
        ${OpenBLAS_SOURCE_DIR}/lapack-netlib/LAPACKE/include)
    endif()
    message("Setting include directory for openblas to ${OpenBLAS_INCLUDE_DIRS}")
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
  SET(LAGRANGE_BUILD_VERSION "v0.5.1-release")
endif()

set(LAGRANGE_DEFINES ${LAGRANGE_DEFINES} LAGRANGE_BUILD_VERSION=${LAGRANGE_BUILD_VERSION})

add_subdirectory(src)

option(BUILD_TESTS "Enable building the tests" OFF)

if(${BUILD_TESTS})
  include(CTest)
  enable_testing()
  add_subdirectory(tests)
endif()
