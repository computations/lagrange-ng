include(FetchContent)
include(ExternalProject)

# nlopt first

#   find_library(nlopt REQUIRED
#           NAMES nlopt libnlopt
#           HINTS "/usr/lib/"
#           )

add_subdirectory(nlopt EXCLUDE_FROM_ALL)

set(BLAS_LIBS_INCLUDE_DIRS ${BLAS_LIBS_INSTALL_PREFIX}/include)
set(BLAS_LIBS_BINARY ${BLAS_LIBS_INSTALL_PREFIX}/lib/libblis.a)
set(BLAS_LIBS_BINARY_SHARED ${BLAS_LIBS_INSTALL_PREFIX}/lib/libblis.so)

set(BLAS_LIBS_INCLUDE_DIRS ${BLAS_LIBS_INCLUDE_DIRS} PARENT_SCOPE)
set(BLAS_LIBS_BINARY ${BLAS_LIBS_BINARY} PARENT_SCOPE)

set(LAPACK OFF CACHE BOOL "Building lapack")
set(LAPACKE ON CACHE BOOL "Building lapacke")
set(USE_OPTIMIZED_BLAS ON CACHE BOOL "Turning off LAPACK blas")

set(BUILD_SINGLE ON CACHE BOOL "Build single precision real" )
set(BUILD_DOUBLE ON CACHE BOOL "Build double precision real" )
set(BUILD_COMPLEX OFF CACHE BOOL "Build single precision complex" )
set(BUILD_COMPLEX16 OFF CACHE BOOL "Build double precision complex" )

set(LAPACKE_BUILD_SINGLE ON CACHE BOOL "Build LAPACKE single precision real" )
set(LAPACKE_BUILD_DOUBLE ON CACHE BOOL "Build LAPACKE double precision real" )
set(LAPACKE_BUILD_COMPLEX OFF CACHE BOOL "Build LAPACKE single precision complex" )
set(LAPACKE_BUILD_COMPLEX16 OFF CACHE BOOL "Build LAPACKE double precision complex" )

set(BLAS_LIBRARIES ${BLAS_LIBS_BINARY_SHARED})
#add_subdirectory(lapack EXCLUDE_FROM_ALL)
add_subdirectory(lapack EXCLUDE_FROM_ALL)

set(LAPACK_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/lapack/LAPACKE/include
  PARENT_SCOPE)

# Jsonhpp fourth

include(FetchContent)

FetchContent_Declare(json
  GIT_REPOSITORY https://github.com/ArthurSonzogni/nlohmann_json_cmake_fetchcontent
  GIT_TAG v3.9.1)

FetchContent_GetProperties(json)
if(NOT json_POPULATED)
  FetchContent_Populate(json)
  add_subdirectory(${json_SOURCE_DIR} ${json_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

#include(cmake/GitUtils.cmake)
#
#git_clone(
#    PROJECT_NAME jsonhpp
#    GIT_URL      https://github.com/nlohmann/json.git
#    GIT_TAG      v3.9.1
#    DIRECTORY    ${CMAKE_CURRENT_SOURCE_DIR}
#    )
#
#add_subdirectory(jsonhpp)
