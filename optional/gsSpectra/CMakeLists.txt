######################################################################
## CMakeLists.txt --- gsSpectra
## This file is part of the G+Smo library.
##
######################################################################

## gsSpectra module
project(gsSpectra)

# Apply G+Smo config
include(gsConfig)

## Collect files
aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR}/src ${PROJECT_NAME}_H)

# Header-only interface library: Spectra is header-only and this module
# adds no .cpp/.hpp of its own, so an OBJECT library (the gsModuleTemplate
# default) would carry an empty $<TARGET_OBJECTS:...>, which is not
# reliably accepted by all generators under this project's CMake floor.
add_library(${PROJECT_NAME} INTERFACE)
target_sources(${PROJECT_NAME} INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_H}>)
target_link_libraries(${PROJECT_NAME} INTERFACE Spectra::Spectra)

#Symlink include dir (in case your headers are in /src)
execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_BINARY_DIR}/${PROJECT_NAME})

# Installed layout mirrors the build-tree symlink, so that
# #include <gsSpectra/gsSpectra.h> resolves in both cases
install(FILES ${${PROJECT_NAME}_H}
  DESTINATION include/gismo/${PROJECT_NAME} )

# Vendored Spectra headers are only installed when this module fetched
# them itself; a system Spectra (found branch) is not ours to redistribute.
if(Spectra_VENDORED)
  get_target_property(_spectra_inc Spectra::Spectra INTERFACE_INCLUDE_DIRECTORIES)
  install(DIRECTORY ${_spectra_inc}/Spectra
          DESTINATION include/gismo
          PATTERN "*.txt" EXCLUDE
          PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)
endif()

# add example files
aux_cpp_directory(${CMAKE_CURRENT_SOURCE_DIR}/examples FILES)
foreach(file ${FILES})
    add_gismo_executable(${file})
    get_filename_component(tarname ${file} NAME_WE) # name without extension
    set_property(TEST ${tarname} PROPERTY LABELS "${PROJECT_NAME}")
    set_target_properties(${tarname} PROPERTIES FOLDER "${PROJECT_NAME}")
    # install the example executables (optionally)
    install(TARGETS ${tarname} DESTINATION "${BIN_INSTALL_DIR}" COMPONENT exe OPTIONAL)
endforeach(file ${FILES})

# add unittests
aux_gs_cpp_directory(${PROJECT_SOURCE_DIR}/unittests unittests_SRCS)
set(gismo_UNITTESTS ${gismo_UNITTESTS} ${unittests_SRCS}
  CACHE INTERNAL "gismo list of unittests")

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)
