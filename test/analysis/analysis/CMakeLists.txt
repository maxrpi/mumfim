add_mpi_test(analysis_volume volume.cc 2 ${MODEL_DIR}/neuron_with_cell.smd ${MESH_DIR}/neuron_with_cell_30k_2_parts.sms 2)
target_link_libraries(analysis_volume util analysis)

add_mpi_test(analysis_measure measure.cc 1)
target_link_libraries(analysis_measure util analysis)
