
add_library(mumfim_torch_utils TorchUtilities.cpp)
set_target_properties(mumfim_torch_utils PROPERTIES PUBLIC_HEADER TorchUtilities.h)
set_target_properties(mumfim_torch_utils PROPERTIES EXPORT_NAME torch_utils)

add_library(mumfim::torch_utils ALIAS mumfim_torch_utils)

find_package(Torch)
find_package(SCOREC REQUIRED)
target_link_libraries(mumfim_torch_utils PUBLIC ${TORCH_LIBRARIES} SCOREC::apf)

target_include_directories(mumfim_torch_utils PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>
        $<INSTALL_INTERFACE:include>)

install(
        TARGETS mumfim_torch_utils
        EXPORT mumfim_torch_utils-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/include)
install(
        EXPORT mumfim_torch_utils-targets
        NAMESPACE mumfim::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/mumfim)