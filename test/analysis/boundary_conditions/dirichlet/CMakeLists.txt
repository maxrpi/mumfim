add_mpi_test(analysis_sim_dirichlet test.cc 1 ${MODEL_DIR}/cube_uniaxial_10kPa.smd ${MESH_DIR}/cube_24_tets.sms)
target_link_libraries(analysis_sim_dirichlet util analysis)

#add_test(NAME dirichlet_bc_regression
#         COMMAND ${CMAKE_COMMAND} -E compare_files "${CMAKE_CURRENT_SOURCE_DIR}/output" "${CMAKE_CURRENT_BINARY_DIR}/output")
