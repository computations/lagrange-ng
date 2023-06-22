cmake_minimum_required(VERSION 3.13)
cmake_policy(SET CMP0135 NEW)

project(lagrange-ng)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

include(FetchContent)

# Find NLopt, otherwise download it
FetchContent_Declare(NLopt
  URL https://github.com/stevengj/nlopt/archive/refs/tags/v2.7.1.tar.gz
  URL_HASH SHA256=db88232fa5cef0ff6e39943fc63ab6074208831dc0031cf1545f6ecd31ae2a1a
  FIND_PACKAGE_ARGS
)

# Find json lib, otherwise download it

FetchContent_Declare(nlohmann_json
  URL https://github.com/ArthurSonzogni/nlohmann_json_cmake_fetchcontent/archive/refs/tags/v3.11.2.tar.gz
  URL_HASH SHA256=8221d9885db8af5c48b472fa2971318443a502c41d8ed6c9bce3a88d6834aa9c
  FIND_PACKAGE_ARGS
)

FetchContent_MakeAvailable(NLopt nlohmann_json)

# Find MKL, fallback to OpenBLAS
set(BLA_VENDOR Intel10_64_dyn)
find_package(BLAS)
find_package(LAPACK)

if (${BLAS_FOUND} AND ${LAPACK_FOUND})
  set(LAGRANGE_LIBS BLAS::BLAS LAPACK::LAPACK)
  set(LAGRANGE_INCLUDES)
  set(LAGRANGE_DEFINES MKL_ENABLED)
else()
    FetchContent_Declare(OpenBLAS
      URL https://github.com/xianyi/OpenBLAS/releases/download/v0.3.23/OpenBLAS-0.3.23.tar.gz
      URL_HASH SHA256=5d9491d07168a5d00116cdc068a40022c3455bf9293c7cb86a65b1054d7e5114
      FIND_PACKAGE_ARGS
    )
  FetchContent_MakeAvailable(OpenBLAS)
    set(LAGRANGE_LIBS blas lapack)
endif()

add_subdirectory(src)

option(BUILD_TESTS "Enable building the tests" OFF)

if(${BUILD_TESTS})
  add_subdirectory(tests)
endif()
