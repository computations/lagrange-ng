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
        IO.cpp
        ConfigFileParser.cpp
        ConfigFileLexer.cpp
        AdjustmentMatrix.cpp
        CSV.cpp
        )

add_library(lagrange-release-libs OBJECT ${LAGRANGE_SOURCES})
add_executable(lagrange-ng main.cpp)

option(LAGRANGE_USER_BARRIER "Force the use of the userspace barrier" OFF)

if(LAGRANGE_USER_BARRIER)
  set(LAGRANGE_DEFINES ${LAGRANGE_DEFINES} LAGRANGE_USER_BARRIER)
endif()

if("${CMAKE_BUILD_TYPE}" MATCHES "Release")
  set(LAGRANGE_BUILD_FLAGS 
    -O3 
    -mavx2 
    -mavx 
    -msse4 
    -mbmi 
    -mbmi2 
    -fno-omit-frame-pointer 
    -g 
    -Wall 
    -Wextra 
    -fdiagnostics-color=always)
    if(${BUILD_ONLY_STATIC})
      target_link_libraries(lagrange-release-libs PUBLIC -static)
    endif()
elseif("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
  set(LAGRANGE_BUILD_FLAGS -g -fno-omit-frame-pointer -fsanitize=address
    -fsanitize-recover=address -fsanitize=undefined -fdiagnostics-color=always)
  set(LAGRANGE_LINK_FLAGS -fsanitize=address -fsanitize=undefined
    -fdiagnostics-color=always)
endif()

set(LAGRANGE_LIBS ${LAGRANGE_LIBS} nlohmann_json::nlohmann_json ${NLOPT_LIB}
  logger)

set(WARNING_FLAGS -Wall -Wextra)

target_compile_options(lagrange-release-libs PUBLIC ${LAGRANGE_BUILD_FLAGS} ${WARNING_FLAGS})
target_include_directories(lagrange-release-libs PUBLIC ${LAGRANGE_INCLUDES})
target_link_libraries(lagrange-release-libs PUBLIC ${LAGRANGE_LIBS} pthread)
target_compile_definitions(lagrange-release-libs PUBLIC ${LAGRANGE_DEFINES})
target_link_options(lagrange-release-libs PUBLIC ${LAGRANGE_LINK_FLAGS})

target_compile_options(lagrange-ng PRIVATE ${LAGRANGE_BUILD_FLAGS} ${WARNING_FLAGS})
target_include_directories(lagrange-ng PRIVATE ${LAGRANGE_INCLUDES})
target_link_libraries(lagrange-ng PRIVATE lagrange-release-libs pthread)
