#add_mpi_test(analysis_error_estimator error_estimator.cc 1 ${MODEL_DIR}/cube_uniaxial_10kPa.smd ${MESH_DIR}/cube_24_tets.sms)
#target_link_libraries(analysis_error_estimator analysis)

#file(COPY        ${CMAKE_CURRENT_SOURCE_DIR}/petsc_options
#     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
