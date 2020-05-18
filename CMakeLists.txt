cmake_minimum_required(VERSION 3.13)
project(lagrange_cpp)
set(CMAKE_CXX_STANDARD 14)

add_subdirectory(deps)
message(STATUS "${ARMADILLO_INCLUDES}")
set(INCLUDE_DIRS /usr/local/include src ${ARMADILLO_INCLUDES})

link_directories(/usr/local/lib)

enable_language(Fortran)

add_executable(lagrange_cpp
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
        src/OptimizeBioGeoAllDispersal.cpp
        src/OptimizeBioGeoAllDispersal.h
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
        src/tree_utils.cpp
        src/tree_utils.h
        src/Utils.cpp
        src/Utils.h
        src/vector_node_object.h
        src/blas.f
        src/clock.f
        src/lapack.f
        src/mataid.f
        src/my_expokit.f
        src/my_matexp.f)

#set(CMAKE_Fortran_FLAGS_DEBUG "-std=legacy")
#set(CMAKE_Fortran_FLAGS "-std=legacy")

set(CMAKE_Fortran_FLAGS_DEBUG "-fallow-argument-mismatch")
set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch")

target_include_directories(lagrange_cpp PRIVATE ${ARMADILLO_INCLUDES})
target_link_libraries(lagrange_cpp gsl nlopt gslcblas ${ARMADILLO_LIBS})
add_dependencies(lagrange_cpp armadillo_headers)
