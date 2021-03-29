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

set(OpenBLAS_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/OpenBLAS)
set(OpenBLAS_INCLUDE_DIRS ${OpenBLAS_INSTALL_PREFIX}/include PARENT_SCOPE)
set(OpenBLAS_LIB ${OpenBLAS_INSTALL_PREFIX}/lib/libopenblas.a PARENT_SCOPE)

include(ExternalProject)
ExternalProject_Add(OpenBLAS
    GIT_REPOSITORY git@github.com:xianyi/OpenBLAS.git
    GIT_TAG v0.3.14
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make
    BUILD_IN_SOURCE TRUE
    INSTALL_COMMAND make PREFIX=${OpenBLAS_INSTALL_PREFIX} install
)
