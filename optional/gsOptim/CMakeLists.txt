######################################################################
## CMakeLists.txt --- gsOptim
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
######################################################################

project(gsOptim)

include( gsConfig)

# aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HEADERS)
# aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_SOURCES)
# aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HPPFILES)
# message ("*** ${PROJECT_NAME} sources: ${${PROJECT_NAME}_HEADERS}, ${${PROJECT_NAME}_SOURCES}, ${${PROJECT_NAME}_HPPFILES}")

#fetch here!
include(gsFetch)
gismo_fetch_directory(optim
  URL https://github.com/kthohr/optim/archive/refs/heads/master.zip
  DESTINATION    external)

IF(NOT EXISTS ${gismo_externals}/optim/include/BaseMatrixOps)
  gismo_fetch_directory(BaseMatrixOps
    URL https://github.com/kthohr/BaseMatrixOps/archive/refs/heads/main.zip
    DESTINATION    external)

  file(RENAME ${gismo_externals}/BaseMatrixOps ${gismo_externals}/optim/include/BaseMatrixOps)
ENDIF()

## CALL ./configure --header-only-version
execute_process (
    COMMAND ./configure --header-only-version
    WORKING_DIRECTORY ${gismo_externals}/optim
)


## CALL ./configure --header-only-version




# Set optim include directory
set(OPTIM_INCLUDE_DIR ${gismo_SOURCE_DIR}/external/optim/header_only_version CACHE INTERNAL "")
# Add optim include directory to G+Smo standard include directories
set(GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${OPTIM_INCLUDE_DIR}
    CACHE INTERNAL "Gismo include directories" FORCE)

message(${GISMO_INCLUDE_DIRS})

include_directories(${GISMO_INCLUDE_DIRS})

add_library(${PROJECT_NAME} INTERFACE)
target_sources(${PROJECT_NAME} INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_HEADERS}>)
add_dependencies(${PROJECT_NAME} optim)

# Set standard properties for all G+Smo extensions
set_target_properties(${PROJECT_NAME} PROPERTIES
  COMPILE_DEFINITIONS gismo_EXPORTS
  POSITION_INDEPENDENT_CODE ON
  LINKER_LANGUAGE CXX
  #START Export all symbols from this extension
  CXX_VISIBILITY_PRESET default
  C_VISIBILITY_PRESET default
  VISIBILITY_INLINES_HIDDEN 0
  #END Export all symbols from this extension
  FOLDER "G+Smo extensions"
  )

# For gsXmlUtils.h
install(DIRECTORY ${CMAKE_SOURCE_DIR}/external/optim/header_only_version
  DESTINATION include/gismo
  PATTERN "*/.svn" EXCLUDE
  )

install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo
        FILES_MATCHING PATTERN "gsOptim/*.h")

# # Add gsXBraid extension to the list of G+Smo extensions
# set(gismo_EXTENSIONS ${gismo_EXTENSIONS} $<TARGET_OBJECTS:${PROJECT_NAME}>
#   CACHE INTERNAL "gismo extensions to be included")

# add example files
#include_directories(${OPTIM_INCLUDE_DIR})
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




# set(${PROJECT_NAME}_HEADERS
# ${gismo_externals}/optim/header_only_version/optim.hpp
# ${gismo_externals}/optim/header_only_version/BaseMatrixOps/include/BaseMatrixOps.hpp
# ${gismo_externals}/optim/header_only_version/constrained/sumt.hpp
# ${gismo_externals}/optim/header_only_version/line_search/more_thuente.hpp
# ${gismo_externals}/optim/header_only_version/misc/bounds_check.hpp
# ${gismo_externals}/optim/header_only_version/misc/determine_bounds_type.hpp
# ${gismo_externals}/optim/header_only_version/misc/error_reporting.hpp
# ${gismo_externals}/optim/header_only_version/misc/error_reporting.ipp
# ${gismo_externals}/optim/header_only_version/misc/jacobian_adjust.hpp
# ${gismo_externals}/optim/header_only_version/misc/numerical_gradient.hpp
# ${gismo_externals}/optim/header_only_version/misc/numerical_hessian.hpp
# ${gismo_externals}/optim/header_only_version/misc/optim_misc.hpp
# ${gismo_externals}/optim/header_only_version/misc/optim_options.hpp
# ${gismo_externals}/optim/header_only_version/misc/optim_structs.hpp
# ${gismo_externals}/optim/header_only_version/misc/optim_trace.hpp
# ${gismo_externals}/optim/header_only_version/misc/transform_vals.hpp
# ${gismo_externals}/optim/header_only_version/stats/optim_stats.hpp
# ${gismo_externals}/optim/header_only_version/stats/seed_values.hpp
# ${gismo_externals}/optim/header_only_version/unconstrained/bfgs.hpp
# ${gismo_externals}/optim/header_only_version/unconstrained/cg.hpp
# ${gismo_externals}/optim/header_only_version/unconstrained/de.hpp
# ${gismo_externals}/optim/header_only_version/unconstrained/de_prmm.hpp
# ${gismo_externals}/optim/header_only_version/unconstrained/gd.hpp
# ${gismo_externals}/optim/header_only_version/unconstrained/lbfgs.hpp
# ${gismo_externals}/optim/header_only_version/unconstrained/newton.hpp
# ${gismo_externals}/optim/header_only_version/unconstrained/nm.hpp
# ${gismo_externals}/optim/header_only_version/unconstrained/optim_unconstrained.hpp
# ${gismo_externals}/optim/header_only_version/unconstrained/pso.hpp
# ${gismo_externals}/optim/header_only_version/unconstrained/pso_dv.hpp
# ${gismo_externals}/optim/header_only_version/zeros/broyden.hpp
# ${gismo_externals}/optim/header_only_version/zeros/broyden_df.hpp
# ${gismo_externals}/optim/header_only_version/zeros/optim_zeros.hpp
# )
