#add_mpi_test(analysis_nonlinear_elastic main.cc 1 ${MODEL_DIR}/cube_uniaxial_10kPa.smd ${MESH_DIR}/cube_24_tets.sms)
#target_link_libraries(analysis_nonlinear_elastic analysis)

#add_test(NAME nonlinear_elastic_regression
#         COMMAND ${CMAKE_COMMAND} -E compare_files "${CMAKE_CURRENT_SOURCE_DIR}/output" "${CMAKE_CURRENT_BINARY_DIR}/output")

#file(COPY        ${CMAKE_CURRENT_SOURCE_DIR}/petsc_options
#     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
