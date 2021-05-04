# Download and unpack googletest at configure time
cmake_minimum_required(VERSION 3.13)
configure_file(deps/CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories("${gtest_SOURCE_DIR}/include")
endif()

set(TREEPATH ${CMAKE_CURRENT_SOURCE_DIR}/data)

# Now simply link against gtest or gtest_main as needed. Eg
set(LAGRANGE_TEST_SOURCES 
  src/tree_reader.cpp
  src/workspace.cpp
  src/operation.cpp
  src/context.cpp
  src/tree.cpp
  src/main.cpp)
add_executable(lagrange-test EXCLUDE_FROM_ALL ${LAGRANGE_TEST_SOURCES} )
target_link_libraries(lagrange-test gtest_main lagrange-release-libs
  ${BLAS_LIBS_BINARY})
target_include_directories(lagrange-test PRIVATE ${CMAKE_SOURCE_DIR}/src/
  ${BLAS_LIBS_INCLUDE_DIRS} ${LAPACK_INCLUDE_DIRS})
target_compile_options(lagrange-test PRIVATE ${LAGRANGE_DEBUG_FLAGS} -fno-omit-frame-pointer -fsanitize=address)
target_link_options(lagrange-test PRIVATE -fno-omit-frame-pointer -fsanitize=address)
target_compile_definitions(lagrange-test PRIVATE -DTREEPATH=${TREEPATH})
