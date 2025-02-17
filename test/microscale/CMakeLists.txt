if(NOT ${BUILD_EXTERNAL})
find_package(amsi REQUIRED)
endif()

add_mpi_test(rve test_rve.cc 1)
target_link_libraries(rve PUBLIC micro_fo)

add_mpi_test(fiber_only ../../src/fiber_only.cc 1 ${CMAKE_CURRENT_BINARY_DIR}/fiber_only.yaml)

add_mpi_test(batched_fiber_only ../../src/batched_fiber_only.cc 1 ${CMAKE_CURRENT_BINARY_DIR}/fiber_only.yaml 1)

file(COPY "fiber_only.yaml" DESTINATION
          "${CMAKE_CURRENT_BINARY_DIR}")

add_mpi_test(test_global_stiffness test_global_stiffness.cc 1)
target_link_libraries(test_global_stiffness PUBLIC micro_fo)

file(COPY "test_global_stiffness_data" DESTINATION
          "${CMAKE_CURRENT_BINARY_DIR}/")

# FIXME This test doesn't actually check anything...need to update this test case to make sure
# that we actually test that the elemental stiffness matrix is getting computed properly!
add_mpi_test(test_elemental_stiffness test_elemental_stiffness.cc 1 ${CMAKE_CURRENT_BINARY_DIR}/fiber_only.yaml)
target_link_libraries(test_elemental_stiffness PUBLIC micro_fo)

#add_mpi_test(test_mass_integrator ./test_mass_integrator.cc 1)
#target_link_libraries(test_mass_integrator micro_fo)

file(COPY fiber_networks
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)


# new unit tests require catch 2
if (NOT ${BUILD_EXTERNAL})
    find_package(Catch2)
endif ()
if (Catch2_FOUND)
    include(Catch)

    set(TEST_SOURCES test_main.cc
            test_reorder_mesh.cc
            test_kokkos_view_of_views.cc
            test_polar_decomposition.cc
            test_finite_difference_material_stiffness.cc
            )
    add_executable(unit_tests ${TEST_SOURCES})
    target_link_libraries(unit_tests PUBLIC Catch2::Catch2 micro_fo)
    catch_discover_tests(unit_tests)
endif ()
