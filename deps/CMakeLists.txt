include(FetchContent)
include(ExternalProject)

# nlopt first

#   find_library(nlopt REQUIRED
#           NAMES nlopt libnlopt
#           HINTS "/usr/lib/"
#           )

add_subdirectory(nlopt EXCLUDE_FROM_ALL)

set(OpenBLAS_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/OpenBLAS)
set(OpenBLAS_INCLUDE_DIRS ${OpenBLAS_INSTALL_PREFIX}/include PARENT_SCOPE)
set(OpenBLAS_LIB_BINARY ${OpenBLAS_INSTALL_PREFIX}/lib/libopenblas.a PARENT_SCOPE)

ExternalProject_Add(OpenBLAS
  GIT_REPOSITORY https://github.com/xianyi/OpenBLAS
  GIT_TAG v0.3.15
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
