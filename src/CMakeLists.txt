add_subdirectory(mumfim/analysis)
add_subdirectory(mumfim/macroscale)
add_subdirectory(mumfim/microscale)

if (NOT ${BUILD_EXTERNAL})
    find_package(amsi REQUIRED)
    find_package(ModelTraits)
    find_package(las)
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
endif ()

add_executable(mumfim_thermal_evaluation thermal_evaluation.cc)
target_link_libraries(mumfim_thermal_evaluation biotissue amsi::amsi LAPACK::LAPACK BLAS::BLAS ModelTraits::ModelTraits ModelTraits::IO MPI::MPI_CXX )
target_compile_features(mumfim_thermal_evaluation PUBLIC cxx_std_17)

add_executable(mumfim_singlescale singlescale_analysis.cc)
target_link_libraries(mumfim_singlescale biotissue amsi::amsi ModelTraits::ModelTraits ModelTraits::IO MPI::MPI_CXX)
target_compile_features(mumfim_singlescale PUBLIC cxx_std_17)
#target_compile_options(mumfim_singlescale PUBLIC $<$<CONFIG:DEBUG>: -fsanitize=address -fno-omit-frame-pointer>)
#target_link_options(mumfim_singlescale PUBLIC $<$<CONFIG:DEBUG>: -fsanitize=address -fno-omit-frame-pointer>)


add_executable(mumfim_multiscale multiscale_analysis.cc)
target_compile_features(mumfim_multiscale PUBLIC cxx_std_17)
target_link_libraries(mumfim_multiscale biotissue micro_fo amsi::amsi ModelTraits::ModelTraits ModelTraits::IO MPI::MPI_CXX)
#target_compile_options(mumfim_multiscale PUBLIC $<$<CONFIG:DEBUG>: -fsanitize=address -fno-omit-frame-pointer>)
#target_link_options(mumfim_multiscale PUBLIC $<$<CONFIG:DEBUG>: -fsanitize=address -fno-omit-frame-pointer>)

add_executable(fiber_only fiber_only.cc)
target_link_libraries(fiber_only PUBLIC micro_fo)

add_executable(batched_fiber_only batched_fiber_only.cc)
target_link_libraries(batched_fiber_only PUBLIC micro_fo)

SET(MUMFIM_INSTALL_TARGETS mumfim_singlescale mumfim_multiscale fiber_only batched_fiber_only)

if (MUMFIM_ENABLE_ADIOS2)
  find_package(adios2 REQUIRED)
  add_executable(adios_fiber_only adios_fiber_only.cc)
  target_link_libraries(adios_fiber_only PUBLIC micro_fo adios2::cxx11_mpi adios2::core adios2::core_mpi adios2::adios2)
  LIST(APPEND MUMFIM_INSTALL_TARGETS adios_fiber_only)
endif()





install(TARGETS ${MUMFIM_INSTALL_TARGETS} RUNTIME DESTINATION bin)
