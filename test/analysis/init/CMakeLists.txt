add_mpi_test(analysis_init test.cc 1)
target_link_libraries(analysis_init analysis)

file(COPY        ${CMAKE_CURRENT_SOURCE_DIR}/petsc_options
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)