# nlopt first

#   find_library(nlopt REQUIRED
#           NAMES nlopt libnlopt
#           HINTS "/usr/lib/"
#           )

include(FetchContent)
include(ExternalProject)

add_subdirectory(nlopt EXCLUDE_FROM_ALL)

# blaze second

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/cmake/blaze_cmake.txt
        DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/blaze)
file(RENAME ${CMAKE_CURRENT_SOURCE_DIR}/blaze/blaze_cmake.txt
        ${CMAKE_CURRENT_SOURCE_DIR}/blaze/CMakeLists.txt)
add_subdirectory(blaze)
set(BLAZE_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/blaze" PARENT_SCOPE)
Blaze_Import(
   QUIET
   BLAS on
   LAPACK on
   THREADING on
   CACHE_SIZE auto
   VECTORIZATION on
   STORAGE_ORDER rowMajor
)

# OpenBLAS third

#find_package(BLAS)
#find_package(LAPACK)
#if(${LAPACK_FOUND} AND ${BLAS_FOUND})
#    set(OpenBLAS_INCLUDE_DIRS lapack blas PARENT_SCOPE)
#    set(OpenBLAS_LIB lapack blas PARENT_SCOPE)
#else()
  #    set(OpenBLAS_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/OpenBLAS)
  #    set(OpenBLAS_INCLUDE_DIRS ${OpenBLAS_INSTALL_PREFIX}/include PARENT_SCOPE)
  #    set(OpenBLAS_LIB ${OpenBLAS_INSTALL_PREFIX}/lib/libopenblas.a PARENT_SCOPE)
  #
  #    include(ExternalProject)
  #    ExternalProject_Add(OpenBLAS
  #        GIT_REPOSITORY https://github.com/xianyi/OpenBLAS
  #        GIT_TAG v0.3.14
  #        GIT_SHALLOW TRUE
  #        GIT_PROGRESS TRUE
  #        CONFIGURE_COMMAND ""
  #        BUILD_COMMAND make USE_OPENMP=1 USE_THREAD=1
  #        BUILD_IN_SOURCE TRUE
  #        INSTALL_COMMAND make PREFIX=${OpenBLAS_INSTALL_PREFIX} install
  #    )
  #    add_dependencies(lagrange-release-libs OpenBLAS)
  #    add_dependencies(lagrange OpenBLAS)
  set(BLAS_LIBS_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/blas)
  set(BLAS_LIBS_INCLUDE_DIRS ${BLAS_LIBS_INSTALL_PREFIX}/include PARENT_SCOPE)
  set(BLAS_LIBS_BINARY ${BLAS_LIBS_INSTALL_PREFIX}/lib/libblis.a PARENT_SCOPE)
  ExternalProject_Add(blis
      GIT_REPOSITORY https://github.com/flame/blis
      GIT_PROGRESS TRUE
      GIT_TAG 0.6.0
      CONFIGURE_COMMAND ./configure --prefix=${BLAS_LIBS_INSTALL_PREFIX} --disable-shared --enable-threading=openmp haswell
      BUILD_COMMAND make
      INSTALL_COMMAND make install
      BUILD_IN_SOURCE TRUE
      UPDATE_DISCONNECTED 1
  )
  add_dependencies(lagrange-release-libs blis)
  add_dependencies(lagrange blis)
  set(LAPACK ON CACHE BOOL "Building lapack")
  set(LAPACKE ON CACHE BOOL "Building lapack")
  set(BUILD_SINGLE OFF CACHE BOOL "Turning off single precision")
  set(BUILD_COMPLEX OFF CACHE BOOL "Turning off single precision complex")
  set(BUILD_COMPLEX16 OFF CACHE BOOL "Turning off double precision complex")
  set(USE_OPTIMIZED_BLAS ON CACHE BOOL "Turning off LAPACK blas")
  add_subdirectory(lapack EXCLUDE_FROM_ALL)
  #endif()

# Jsonhpp fourth


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
