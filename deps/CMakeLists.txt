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
