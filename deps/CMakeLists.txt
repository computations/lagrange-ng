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
