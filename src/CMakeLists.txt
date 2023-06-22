set(LAGRANGE_SOURCES
        AncSplit.cpp
        Context.cpp
        Node.cpp
        Tree.cpp
        TreeReader.cpp
        Utils.cpp
        Workspace.cpp
        Operation.cpp
        WorkerState.cpp
        Arnoldi.cpp
        Fossil.cpp
        ConfigFile.cpp
        Alignment.cpp
        Periods.cpp
        )

add_library(lagrange-release-libs OBJECT ${LAGRANGE_SOURCES})
add_executable(lagrange main.cpp)

if("${CMAKE_BUILD_TYPE}" MATCHES "Release")
  set(LAGRANGE_BUILD_FLAGS -O3 -mtune=native -march=native -fno-omit-frame-pointer)
elseif("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
  set(LAGRANGE_BUILD_FLAGS -g -fno-omit-frame-pointer -fsanitize=address -fsanitize-recover=address)
  set(LAGRANGE_LINK_FLAGS -fsanitize=address)
endif()


set(LAGRANGE_LIBS nlohmann_json::nlohmann_json NLopt::nlopt ${MKL_LIBRARIES})

target_compile_options(lagrange-release-libs PUBLIC ${LAGRANGE_BUILD_FLAGS} ${WARNING_FLAGS})
target_include_directories(lagrange-release-libs PUBLIC ${MKL_INCLUDE_DIRS} ${MKL_INCLUDE_DIRS})
target_link_libraries(lagrange-release-libs PUBLIC ${LAGRANGE_LIBS} pthread)
target_compile_definitions(lagrange-release-libs PUBLIC MKL_ENABLED)
# target_link_directories(lagrange-release-libs PUBLIC ${LAGRANGE_LINALG_LINK_DIRS})
# target_link_options(lagrange-release-libs PUBLIC ${LAGRANGE_LINK_FLAGS})

target_compile_options(lagrange PRIVATE ${LAGRANGE_BUILD_FLAGS} ${WARNING_FLAGS})
target_include_directories(lagrange PRIVATE ${MKL_INCLUDE_DIRS})
target_link_libraries(lagrange PRIVATE lagrange-release-libs pthread)
target_compile_definitions(lagrange PRIVATE MKL_ENABLED)
#target_link_directories(lagrange PRIVATE ${LAGRANGE_LINALG_LINK_DIRS})
#target_link_options(lagrange PRIVATE ${LAGRANGE_LINK_FLAGS})
