add_mpi_test(analysis_loading test.cc 1 ${MODEL_DIR}/cube_uniaxial_10kPa.smd ${MESH_DIR}/cube_24_tets.sms)
target_link_libraries(analysis_loading analysis)