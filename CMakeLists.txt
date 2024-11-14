cmake_minimum_required(VERSION 3.13)
cmake_policy(SET CMP0135 NEW)


project(lagrange-ng)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")


include(FetchContent)

# Find NLopt, otherwise download it
FetchContent_Declare(NLopt
  URL https://github.com/stevengj/nlopt/archive/refs/tags/v2.7.1.tar.gz
  URL_HASH SHA256=db88232fa5cef0ff6e39943fc63ab6074208831dc0031cf1545f6ecd31ae2a1a
  EXCLUDE_FROM_ALL
  FIND_PACKAGE_ARGS
)


# Find json lib, otherwise download it
FetchContent_Declare(nlohmann_json
  URL https://github.com/ArthurSonzogni/nlohmann_json_cmake_fetchcontent/archive/refs/tags/v3.11.2.tar.gz
  URL_HASH SHA256=8221d9885db8af5c48b472fa2971318443a502c41d8ed6c9bce3a88d6834aa9c
  EXCLUDE_FROM_ALL
  FIND_PACKAGE_ARGS
)

FetchContent_Declare(logger
  URL https://github.com/computations/logger/archive/refs/tags/v0.2.0.tar.gz
  URL_HASH SHA256=92adc9a8f96c8077f09ac4a97c8c0d28b990579b95b1aad0c16724d4c7fb95bd
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
      URL https://github.com/xianyi/OpenBLAS/releases/download/v0.3.23/OpenBLAS-0.3.23.tar.gz
      URL_HASH SHA256=5d9491d07168a5d00116cdc068a40022c3455bf9293c7cb86a65b1054d7e5114
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

add_subdirectory(src)

option(BUILD_TESTS "Enable building the tests" OFF)

if(${BUILD_TESTS})
  include(CTest)
  enable_testing()
  add_subdirectory(tests)
endif()
