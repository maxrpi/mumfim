find_package(PkgConfig REQUIRED)
pkg_check_modules(PETSc REQUIRED IMPORTED_TARGET PETSc)

find_package(SCOREC CONFIG REQUIRED)

if (NOT ${BUILD_EXTERNAL})
    find_package(amsi REQUIRED)
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/amsiAnalysisConfig.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/amsiAnalysisConfig.h
        @ONLY)

set(analysis_src_files
        amsiAnalysis.cc
        amsiBoundaryConditions.cc
        amsiBoundaryConditionQuery.cc
        amsiDeformation.cc
        amsiFields.cc
        amsiFEA.cc
        amsiLinearElasticConstitutive.cc
        amsiNonlinearAnalysis.cc
        amsiPETScLAS.cc
        apfElementalSystem.cc
        apfFunctions.cc
        apfMatrixUtil.cc
        apfMeasure.cc
        apfMeshIterator.cc
        apfMeshUtil.cc
        apfLabelRegions.cc
        Solvers.cc
        ElementalSystem.cc
        amsiNeumannIntegratorsMT.cc)

# Add the library to create
add_library(analysis ${analysis_src_files})
add_library(mumfim::analysis ALIAS analysis)

# properties for downstream packages
target_include_directories(analysis PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)

if (NOT ${BUILD_EXTERNAL})
    find_package(ModelTraits REQUIRED)
endif ()

target_compile_features(analysis PUBLIC cxx_std_17)

target_link_libraries(analysis PUBLIC amsi::util SCOREC::gmi SCOREC::apf SCOREC::ma PkgConfig::PETSc ModelTraits::ModelTraits)
