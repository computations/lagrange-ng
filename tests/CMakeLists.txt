cmake_minimum_required(VERSION 3.13)
include(FetchContent)

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.17.0.tar.gz
  DOWNLOAD_EXTRACT_TIMESTAMP NEW
)

FetchContent_MakeAvailable(googletest)

set(TREEPATH ${CMAKE_CURRENT_SOURCE_DIR}/data)

# Now simply link against gtest or gtest_main as needed. Eg
set(LAGRANGE_TEST_SOURCES 
  src/tree_reader.cpp
  src/workspace.cpp
  src/operation.cpp
  src/context.cpp
  src/tree.cpp
  src/configfile.cpp
  src/csv.cpp
  src/adjustment_matrix.cpp
  src/main.cpp
)

add_executable(lagrange-test ${LAGRANGE_TEST_SOURCES} )
target_link_libraries(lagrange-test GTest::gtest_main lagrange-release-libs
  ${BLAS_LIBS_BINARY})
target_include_directories(lagrange-test PRIVATE ${CMAKE_SOURCE_DIR}/src/
  ${BLAS_LIBS_INCLUDE_DIRS} ${LAPACK_INCLUDE_DIRS})
target_compile_options(lagrange-test PRIVATE ${LAGRANGE_DEBUG_FLAGS} -fno-omit-frame-pointer -fsanitize=address -fsanitize-recover=address -g -O0)
target_link_options(lagrange-test PRIVATE -fno-omit-frame-pointer -fsanitize=address -fsanitize-recover=address -g -O0)
target_compile_definitions(lagrange-test PRIVATE -DTREEPATH=${TREEPATH})

include(GoogleTest)
gtest_discover_tests(lagrange-test)
