cmake_minimum_required(VERSION 3.13)
project(lagrange-cpp)
set(CMAKE_CXX_STANDARD 14)

set(INCLUDE_DIRS /usr/local/include src ${ARMADILLO_INCLUDES})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

link_directories(/usr/local/lib)
find_package(OpenMP)

set(LAGRANGE_SOURCES
        src/AncSplit.cpp
        src/Context.cpp
        src/InputReader.cpp
        src/Node.cpp
        src/Tree.cpp
        src/TreeReader.cpp
        src/Utils.cpp
        src/Workspace.cpp
        src/Operation.cpp
        src/WorkerState.cpp
        )

set(WARNING_FLAGS -Wall -Wextra -Wshadow)

set(CMAKE_CXX_FLAGS_RELEASE "-O2 -mtune=native -march=native -fno-omit-frame-pointer")
set(CMAKE_C_FLAGS_RELEASE "-O2 -mtune=native -march=native -fno-omit-frame-pointer")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -fno-omit-frame-pointer -fsanitize=address -fsanitize-recover=address -fsanitize=thread")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g -fno-omit-frame-pointer -fsanitize=address -fsanitize-recover=address -fsanitize=thread")
set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address -fsanitize-recover=address -fsanitize=thread")

set(OPT_LIBS nlopt)
set(LAGRANGE_LIBS ${OPT_LIBS} nlohmann_json::nlohmann_json)

add_library(lagrange-release-libs OBJECT ${LAGRANGE_SOURCES})
add_executable(lagrange src/main.cpp)

add_subdirectory(deps)

target_compile_options(lagrange-release-libs PRIVATE ${LAGRANGE_RELEASE_FLAGS} ${WARNING_FLAGS})
target_include_directories(lagrange-release-libs PRIVATE ${OpenBLAS_INCLUDE_DIRS})
target_link_libraries(lagrange-release-libs PRIVATE ${LAGRANGE_LIBS} pthread OpenMP::OpenMP_CXX ${OpenBLAS_LIB_BINARY})

target_compile_options(lagrange PRIVATE ${LAGRANGE_RELEASE_FLAGS} ${WARNING_FLAGS})
target_include_directories(lagrange PRIVATE ${OpenBLAS_INCLUDE_DIRS})
target_link_libraries(lagrange PRIVATE lagrange-release-libs pthread nlohmann_json::nlohmann_json
    OpenMP::OpenMP_CXX ${OpenBLAS_LIB_BINARY})

add_subdirectory(tests)
