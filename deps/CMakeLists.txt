include(FetchContent)
include(ExternalProject)

# nlopt first

#   find_library(nlopt REQUIRED
#           NAMES nlopt libnlopt
#           HINTS "/usr/lib/"
#           )

add_subdirectory(nlopt EXCLUDE_FROM_ALL)

# Jsonhpp second

include(FetchContent)

FetchContent_Declare(json
  GIT_REPOSITORY https://github.com/ArthurSonzogni/nlohmann_json_cmake_fetchcontent
  GIT_TAG v3.9.1)

FetchContent_GetProperties(json)
if(NOT json_POPULATED)
  FetchContent_Populate(json)
  add_subdirectory(${json_SOURCE_DIR} ${json_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()
