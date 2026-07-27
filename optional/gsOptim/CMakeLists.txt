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
  GIT_REPOSITORY https://github.com/hverhelst/optim.git
  DESTINATION    external)

IF(NOT EXISTS ${gismo_externals}/optim/include/BaseMatrixOps/include)
    execute_process(COMMAND "${GIT_EXECUTABLE}" "submodule" "update" "--init"
	                   WORKING_DIRECTORY ${gismo_SOURCE_DIR}/external/optim)
ENDIF()

execute_process (
    #COMMAND dos2unix ./configure
    COMMAND ./configure --header-only-version
    WORKING_DIRECTORY ${gismo_externals}/optim
    RESULT_VARIABLE ret
    #in case of error executing due to windows file: dos2unix script
    )

if(ret AND NOT ret EQUAL 0)
    message( FATAL_ERROR "Configure of optim returned an error.\n(try dos2unix external/optim/configure if on linux)")
  else()
    #BaseMatrixOps is an external dependency of optim
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${gismo_SOURCE_DIR}/external/optim/header_only_version/BaseMatrixOps ${gismo_SOURCE_DIR}/external/BaseMatrixOps)
  endif()

# Set optim include directory
set(OPTIM_INCLUDE_DIR ${gismo_SOURCE_DIR}/external/optim/header_only_version CACHE INTERNAL "")
# Add optim include directory to G+Smo standard include directories
set(GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${OPTIM_INCLUDE_DIR}
    CACHE INTERNAL "Gismo include directories" FORCE)

include_directories(${GISMO_INCLUDE_DIRS})

if(GISMO_BUILD_LIB AND GISMO_BUILD_MODULE_LIB)

  # Runtime-loadable module library: gsOptim_module.so with the ABI
  # handshake entry points; loaded via gsModuleManager
  include(gsModuleLib)
  gismo_add_module_lib(${PROJECT_NAME} "1.0.0"
    ${CMAKE_CURRENT_SOURCE_DIR}/gsOptimModuleEntry.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/gsOptimRegister.cpp)

else()

  # Legacy path: compiled into libgismo; the static registrar in
  # gsOptimRegister.cpp populates the optimizer registry on load
  add_library(${PROJECT_NAME} INTERFACE)
  target_sources(${PROJECT_NAME} INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_HEADERS}>)
  add_dependencies(${PROJECT_NAME} optim)

  set(gismo_EXTENSIONS ${gismo_EXTENSIONS}
      ${CMAKE_CURRENT_SOURCE_DIR}/gsOptimRegister.cpp
      CACHE INTERNAL "Gismo extensions to be included")

endif()

# For gsXmlUtils.h
install(DIRECTORY ${CMAKE_SOURCE_DIR}/external/optim/header_only_version
  DESTINATION include/gismo
  PATTERN "*/.svn" EXCLUDE
  )

install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo
        FILES_MATCHING PATTERN "gsOptim/*.h")

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
