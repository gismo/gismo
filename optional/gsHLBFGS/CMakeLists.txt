######################################################################
### CMakeLists.txt --- gsHLBFGS
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
######################################################################

## gsHLBFGS module
project(gsHLBFGS)

# Apply G+Smo config
include(gsConfig)

## Collect files
aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR}/src ${PROJECT_NAME}_H)

# Interface library: the compiled HLBFGS sources come from
# gismo_add_dependency() in dependencies.cmake, whose objects already
# land in gismo_EXTENSIONS. This module contributes no .cpp of its own,
# so an OBJECT library (the gsModuleTemplate default) would carry an
# empty $<TARGET_OBJECTS:...>, which is not reliably accepted by all
# generators under this project's CMake floor.
add_library(${PROJECT_NAME} INTERFACE)
target_sources(${PROJECT_NAME} INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_H}>)
target_link_libraries(${PROJECT_NAME} INTERFACE HLBFGS::HLBFGS)

#Symlink include dir (in case your headers are in /src)
execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_BINARY_DIR}/${PROJECT_NAME})

# Installed layout mirrors the build-tree symlink, so that
# #include <gsHLBFGS/gsHLBFGS.h> resolves in both cases
install(FILES ${${PROJECT_NAME}_H}
  DESTINATION include/gismo/${PROJECT_NAME} )

# The vendored HLBFGS headers are part of this module's public interface:
# src/gsHLBFGS.h includes "HLBFGS/HLBFGS.h", which resolves against the
# consumer's -I<prefix>/include/gismo, so they install beside gsHLBFGS/,
# not inside it. A non-vendored HLBFGS is not ours to redistribute.
if(HLBFGS_VENDORED)
  get_target_property(_hlbfgs_inc HLBFGS::HLBFGS INTERFACE_INCLUDE_DIRECTORIES)
  install(DIRECTORY ${_hlbfgs_inc}/HLBFGS
          DESTINATION include/gismo
          FILES_MATCHING PATTERN "*.h")
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
