######################################################################
## gsLibrary.cmake
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
######################################################################

#include (GenerateExportHeader)

#message("Using ${${PROJECT_NAME}_EXTENSIONS}")
#message("Using ${${PROJECT_NAME}_MODULES}")
#message("Using ${${PROJECT_NAME}_SOURCES}")

###################################################################
# Static library
###################################################################

add_library(${PROJECT_NAME}_static STATIC
  ${${PROJECT_NAME}_MODULES}
  ${${PROJECT_NAME}_SOURCES}
  ${${PROJECT_NAME}_EXTENSIONS}
  )

#generate_export_header(${PROJECT_NAME}_static)

if(${PROJECT_NAME}_LINKER)
  target_link_libraries(${PROJECT_NAME}_static "${${PROJECT_NAME}_LINKER}")
endif()

if (GISMO_WITH_XDEBUG AND DBGHELP_FOUND)
  target_link_libraries(${PROJECT_NAME}_static ${DBGHELP_LIBRARY})
ENDIF()

if (GISMO_GCC_STATIC_LINKAGE)
  target_link_libraries(${PROJECT_NAME}_static -static-libgcc -static-libstdc++)
endif()

# Avoid naming conflic on MSVC
if("x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xMSVC")
  set(gs_static_lib_suffix _static)
endif()

set_target_properties(${PROJECT_NAME}_static PROPERTIES
COMPILE_DEFINITIONS ${PROJECT_NAME}_STATIC
POSITION_INDEPENDENT_CODE ON
LINKER_LANGUAGE CXX
CXX_VISIBILITY_PRESET "hidden"
FOLDER "G+Smo libraries"
OUTPUT_NAME ${PROJECT_NAME}${gs_static_lib_suffix} )

###################################################################
# Pygismo
###################################################################

if (GISMO_WITH_PYBIND11)

  pybind11_add_module(py${PROJECT_NAME} MODULE
    "${gismo_SOURCE_DIR}/src/misc/gsPyBind11.cpp"
    )

  set_target_properties(py${PROJECT_NAME} PROPERTIES
    COMPILE_DEFINITIONS gismo_EXPORTS
    POSITION_INDEPENDENT_CODE ON
    LINKER_LANGUAGE CXX
    CXX_VISIBILITY_PRESET "hidden"
    CUDA_VISIBILITY_PRESET "hidden"
    )

  # since gismo (${PROJECT_NAME}) target includes bindings, it needs
  # pybind/python info. Those are automatically managed in
  # `pybind11_add_module`. Since we aren't using it, setup gismo target
  # in similar fashion manually. 
  target_link_libraries(${PROJECT_NAME}_static pybind11::module)

  if(NOT DEFINED CMAKE_INTERPROCEDURAL_OPTIMIZATION)
    target_link_libraries(${PROJECT_NAME}_static pybind11::lto)
  endif()

  # link gismo to pygismo
  target_link_libraries(py${PROJECT_NAME} PRIVATE ${PROJECT_NAME}_static)

  pybind11_strip(py${PROJECT_NAME})
  pybind11_extension(py${PROJECT_NAME})

  if (GISMO_KLSHELL)
    target_compile_definitions(py${PROJECT_NAME} PUBLIC GISMO_KLSHELL)
  endif()# To fix

endif(GISMO_WITH_PYBIND11)

###################################################################
# Shared library
###################################################################

if(GISMO_BUILD_LIB)

  #if ("x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xGNU")
  #  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-implicit-templates")
  #endif()

  if("x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xMSVC" OR
     "x${CMAKE_GENERATOR}" STREQUAL "xXcode")
    set(${PROJECT_NAME}_SOURCES ${${PROJECT_NAME}_SOURCES}
      "${gismo_SOURCE_DIR}/src/misc/gsDllMain.cpp")
  endif()

  if (GISMO_WITH_XDEBUG)
    if (NOT "x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xMSVC" OR DBGHELP_FOUND)
      set(${PROJECT_NAME}_SOURCES ${${PROJECT_NAME}_SOURCES} ${gismo_SOURCE_DIR}/src/misc/gsStackWalker.cpp)
    endif()
  endif()

  add_library(${PROJECT_NAME} SHARED
    ${${PROJECT_NAME}_MODULES}
    ${${PROJECT_NAME}_SOURCES}
    ${${PROJECT_NAME}_EXTENSIONS}
    )

  set_target_properties(${PROJECT_NAME} PROPERTIES
    #https://community.kde.org/Policies/Binary_Compatibility_Issues_With_C%2B%2B
    VERSION "${${PROJECT_NAME}_VERSION}"
    SOVERSION "${${PROJECT_NAME}_VERSION_MAJOR}"
    PUBLIC_HEADER "${PROJECT_SOURCE_DIR}/src/${PROJECT_NAME}.h"
    POSITION_INDEPENDENT_CODE ON
    LINKER_LANGUAGE CXX
    CXX_VISIBILITY_PRESET "hidden"
    #COMPILE_DEFINITIONS ${PROJECT_NAME}_EXPORTS # Used for DLL exporting (defined by default by CMake)
    FOLDER "G+Smo libraries"
    )
    #generate_export_header(${PROJECT_NAME})

    if (GISMO_WITH_PYBIND11)
      target_link_libraries(${PROJECT_NAME} pybind11::embed)
    endif()

  #if(gsMpfr_ENABLED OR gsGmp_ENABLED)
  #    find_package(GMP)
  #    find_package(MPFR)
  #
  #    if (GMP_FOUND AND MPFR_FOUND)
  #      target_link_libraries(${PROJECT_NAME} ${MPFR_LIBRARY};${GMP_LIBRARY};${GMPXX_LIBRARY})
  #    endif()
  #endif()

  if (GISMO_WITH_SUPERLU)
    target_link_libraries(${PROJECT_NAME} ${SUPERLU_LIBRARIES})
  endif()

  if (GISMO_WITH_TAUCS)
    target_link_libraries(${PROJECT_NAME} ${TAUCS_LIBRARIES})
  endif()

  if (GISMO_WITH_UMFPACK)
    target_link_libraries(${PROJECT_NAME} ${UMFPACK_LIBRARIES})
  endif()

  if (GISMO_WITH_PARDISO)
     if (PARDISO_USE_MKL)
       find_package(MKL REQUIRED)
       target_link_libraries(${PROJECT_NAME} ${MKL_LIBRARIES})
     else()
       find_package(Pardiso REQUIRED)
       target_link_libraries(${PROJECT_NAME} Pardiso)
     endif()
  endif()

  if(${PROJECT_NAME}_LINKER)
    target_link_libraries(${PROJECT_NAME} "${${PROJECT_NAME}_LINKER}")
  endif()

  if (GISMO_GCC_STATIC_LINKAGE)
    target_link_libraries(${PROJECT_NAME} -static-libgcc -static-libstdc++)
  endif()

  #  if (GISMO_WITH_OPENMP)
  #    find_package(OpenMP REQUIRED)
  #  endif()

  if (GISMO_WITH_XDEBUG AND DBGHELP_FOUND)
    target_link_libraries(${PROJECT_NAME} ${DBGHELP_LIBRARY})
  endif()

  if( WIN32 ) # Copy the dll to the bin folder to allow executables to find it
    if(CMAKE_CONFIGURATION_TYPES)
      add_custom_command(
      TARGET ${PROJECT_NAME}
      POST_BUILD
      #COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/bin
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/bin/$<CONFIGURATION>
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:${PROJECT_NAME}> ${CMAKE_BINARY_DIR}/bin/$<CONFIGURATION>
      COMMAND ${CMAKE_COMMAND} -E echo 'The file $<TARGET_FILE:${PROJECT_NAME}> is copied to the bin folder for convenience.' )
    else()
      add_custom_command(
      TARGET ${PROJECT_NAME}
      POST_BUILD
      #COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/bin
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/bin
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:${PROJECT_NAME}> ${CMAKE_BINARY_DIR}/bin
      COMMAND ${CMAKE_COMMAND} -E echo 'The file $<TARGET_FILE:${PROJECT_NAME}> is copied to the bin folder for convenience.' )
    endif()
  endif( WIN32 )

endif(GISMO_BUILD_LIB)

if (EIGEN_USE_MKL_ALL)
  target_link_libraries(${PROJECT_NAME} ${MKL_LIBRARIES})
endif()

## #################################################################
## Installation
## #################################################################

install(TARGETS ${PROJECT_NAME}_static OPTIONAL
  EXPORT gismoTargets
  LIBRARY DESTINATION "${LIB_INSTALL_DIR}" COMPONENT shared
  ARCHIVE DESTINATION "${LIB_INSTALL_DIR}" COMPONENT static
  RUNTIME DESTINATION "${BIN_INSTALL_DIR}" COMPONENT exe
  PUBLIC_HEADER DESTINATION "${INCLUDE_INSTALL_DIR}/${PROJECT_NAME}")

if(GISMO_BUILD_LIB)

  install(TARGETS ${PROJECT_NAME}
    # IMPORTANT: Add the ${PROJECT_NAME} library to the "export-set"
    EXPORT gismoTargets
    LIBRARY DESTINATION "${LIB_INSTALL_DIR}" COMPONENT shared
    ARCHIVE DESTINATION "${LIB_INSTALL_DIR}" COMPONENT static
    RUNTIME DESTINATION "${BIN_INSTALL_DIR}" COMPONENT exe
    PUBLIC_HEADER DESTINATION "${INCLUDE_INSTALL_DIR}/${PROJECT_NAME}")

endif(GISMO_BUILD_LIB)
