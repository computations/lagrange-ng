# nlopt first

#   find_library(nlopt REQUIRED
#           NAMES nlopt libnlopt
#           HINTS "/usr/lib/"
#           )

add_subdirectory(nlopt)

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
   THREADING off
   CACHE_SIZE auto
   VECTORIZATION on
   STORAGE_ORDER rowMajor
)

# OpenBLAS third
find_package(BLAS)
find_package(LAPACK)
if(${LAPACK_FOUND} AND ${BLAS_FOUND})
    set(OpenBLAS_INCLUDE_DIRS lapack blas PARENT_SCOPE)
    set(OpenBLAS_LIB lapack blas PARENT_SCOPE)
else()
    set(OpenBLAS_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/OpenBLAS)
    set(OpenBLAS_INCLUDE_DIRS ${OpenBLAS_INSTALL_PREFIX}/include PARENT_SCOPE)
    set(OpenBLAS_LIB ${OpenBLAS_INSTALL_PREFIX}/lib/libopenblas.a PARENT_SCOPE)

    include(ExternalProject)
    ExternalProject_Add(OpenBLAS
        GIT_REPOSITORY https://github.com/xianyi/OpenBLAS
        GIT_TAG v0.3.14
        GIT_SHALLOW TRUE
        GIT_PROGRESS TRUE
        CONFIGURE_COMMAND ""
        BUILD_COMMAND make
        BUILD_IN_SOURCE TRUE
        INSTALL_COMMAND make PREFIX=${OpenBLAS_INSTALL_PREFIX} install
    )
    add_dependencies(lagrange-release-libs OpenBLAS)
    add_dependencies(lagrange OpenBLAS)
endif()

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
