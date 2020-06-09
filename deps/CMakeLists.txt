#armadillo first
add_custom_target(armadillo_headers
        COMMAND tar xf
        ${CMAKE_CURRENT_SOURCE_DIR}/armadillo-3.0.2.tar.gz -C
        ${CMAKE_CURRENT_BINARY_DIR})

message(STATUS "${CMAKE_CURRENT_BINARY_DIR}/armadillo-3.0.2/include")

set(ARMADILLO_INCLUDES "${CMAKE_CURRENT_BINARY_DIR}/armadillo-3.0.2/include"
        PARENT_SCOPE)

#nlopt second

find_library(nlopt REQUIRED
        NAMES nlopt libnlopt
        HINTS "/usr/lib/"
        )

# blaze

#add_subdirectory(blaze)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/cmake/blaze_cmake.txt
        DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/blaze)
file(RENAME ${CMAKE_CURRENT_SOURCE_DIR}/blaze/blaze_cmake.txt
        ${CMAKE_CURRENT_SOURCE_DIR}/blaze/CMakeLists.txt)
add_subdirectory(blaze)
set(BLAZE_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/blaze" PARENT_SCOPE)
