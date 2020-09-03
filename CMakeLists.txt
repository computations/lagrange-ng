cmake_minimum_required(VERSION 3.13)
project(lagrange_cpp)
set(CMAKE_CXX_STANDARD 14)

add_subdirectory(deps)
message(STATUS "${ARMADILLO_INCLUDES}")
set(INCLUDE_DIRS /usr/local/include src ${ARMADILLO_INCLUDES})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

link_directories(/usr/local/lib)

enable_language(Fortran)

set(LAGRANGE_SOURCES
        src/AncSplit.cpp
        src/AncSplit.h
        src/BayesianBioGeo.cpp
        src/BayesianBioGeo.h
        src/BayesianBioGeoAllDispersal.cpp
        src/BayesianBioGeoAllDispersal.h
        src/BioGeoTree.cpp
        src/BioGeoTree.h
        src/BioGeoTreeTools.cpp
        src/BioGeoTreeTools.h
        src/BranchSegment.cpp
        src/BranchSegment.h
        src/double_node_object.h
        src/expm.h
        src/InputReader.cpp
        src/InputReader.h
        src/main.cpp
        src/node.cpp
        src/node.h
        src/node_object.h
        src/OptimizeBioGeo.cpp
        src/OptimizeBioGeo.h
        src/OptimizeBioGeoAllDispersal_nlopt.cpp
        src/OptimizeBioGeoAllDispersal_nlopt.h
        src/RateMatrixUtils.cpp
        src/RateMatrixUtils.h
        src/RateModel.cpp
        src/RateModel.h
        src/string_node_object.h
        src/superdouble.cpp
        src/superdouble.h
        src/tree.cpp
        src/tree.h
        src/tree_reader.cpp
        src/tree_reader.h
        src/Utils.cpp
        src/Utils.h
        src/vector_node_object.h
        src/blas.f
        src/clock.f
        src/lapack.f
        src/mataid.f
        src/my_expokit.f
        src/my_matexp.f)


add_executable(lagrange ${LAGRANGE_SOURCES})
add_executable(lagrange_debug ${LAGRANGE_SOURCES})
add_executable(lagrange_profile ${LAGRANGE_SOURCES})

# These flags are required for the expmkit code
set(CMAKE_Fortran_FLAGS_DEBUG "-fallow-argument-mismatch")
set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch")

set(WARNING_FLAGS -Wall -Wextra -Wshadow)

set(LAGRANGE_DEBUG_FLAGS -Og -g -fno-omit-frame-pointer)
list(APPEND ${LAGRANGE_DEBUG_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG})
set(LAGRANGE_PROFILE_FLAGS -O3 -fno-omit-frame-pointer)
list(APPEND ${LAGRANGE_PROFILE_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE})
set(LAGRANGE_RELEASE_FLAGS -O3 -DBLAZE_USE_VECTORIZATION=1 -g)
list(APPEND ${LAGRANGE_RELEASE_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE})

set(LINALG_LIBS gsl cblas lapack)
set(OPT_LIBS nlopt)
set(LAGRANGE_LIBS ${LINALG_LIBS} ${OPT_LIBS})

target_include_directories(lagrange PRIVATE ${BLAZE_INCLUDES})
target_link_libraries(lagrange ${LAGRANGE_LIBS})
target_compile_options(lagrange PRIVATE ${LAGRANGE_RELEASE_FLAGS}
        ${WARNING_FLAGS})

target_include_directories(lagrange_debug PRIVATE ${BLAZE_INCLUDES})
target_link_libraries(lagrange_debug ${LAGRANGE_LIBS})
target_compile_options(lagrange_debug PRIVATE ${LAGRANGE_DEBUG_FLAGS}
        ${WARNING_FLAGS})

target_include_directories(lagrange_profile PRIVATE ${BLAZE_INCLUDES})
target_link_libraries(lagrange_profile ${LAGRANGE_LIBS})
target_compile_options(lagrange_profile PRIVATE ${LAGRANGE_PROFILE_FLAGS}
        ${WARNING_FLAGS})
