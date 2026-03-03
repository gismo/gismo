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
LABELS "${PROJECT_NAME}"
OUTPUT_NAME ${PROJECT_NAME}${gs_static_lib_suffix} )

###################################################################
# Pygismo
###################################################################

if (GISMO_WITH_NANOBIND)

  # Determine the best Python target once
  set(_py_target "")
  if(TARGET Python::Module)
    set(_py_target Python::Module)
  elseif(TARGET Python::Python)
    set(_py_target Python::Python)
  endif()

  # Linker fix for strict CI environments (manlylinux/gcc-toolset)
  # This tells the linker to allow Python symbols to be resolved at runtime
  set(_py_linker_fix "")
  if(NOT MSVC)
    set(_py_linker_fix "-Wl,--allow-shlib-undefined" "-Wl,-z,undefs")
  endif()

  set(PYGISMO_NAME_MAP_gsCore        "core")
  set(PYGISMO_NAME_MAP_gsIO          "io")
  set(PYGISMO_NAME_MAP_gsNurbs       "nurbs")
  set(PYGISMO_NAME_MAP_gsModeling    "modeling")
  set(PYGISMO_NAME_MAP_gsPde         "pde")
  set(PYGISMO_NAME_MAP_gsMatrix      "matrix")
  set(PYGISMO_NAME_MAP_gsHSplines    "hsplines")
  set(PYGISMO_NAME_MAP_gsAssembler   "assembler")
  set(PYGISMO_NAME_MAP_gsMSplines    "msplines")

  nanobind_add_module(pygismo__core
    NB_SHARED
    NB_DOMAIN gismo
    "${gismo_SOURCE_DIR}/src/misc/gsNanoBind.cpp"
  )
  target_link_libraries(pygismo__core PRIVATE ${PROJECT_NAME} ${_py_target})
  target_link_options(pygismo__core PRIVATE ${_py_linker_fix})

  set_target_properties(pygismo__core PROPERTIES
    OUTPUT_NAME "_core"
    LIBRARY_OUTPUT_DIRECTORY "${PYGISMO_PKG_DIR}"
  )

  list(APPEND PYGISMO_TARGETS pygismo__core)

  file(GLOB _nb_binding_files "${gismo_SOURCE_DIR}/src/gs*/nanobind/*_nb.cpp")
  foreach(_nb_file ${_nb_binding_files})
    get_filename_component(_nb_name ${_nb_file} NAME_WE)
    string(REGEX REPLACE "_nb$" "" _src_name "${_nb_name}")

    if(DEFINED PYGISMO_NAME_MAP_${_src_name})
      set(_mod_name "${PYGISMO_NAME_MAP_${_src_name}}")
    else()
      string(TOLOWER "${_src_name}" _mod_name)
    endif()

    nanobind_add_module(pygismo_${_src_name}
      NB_SHARED
      NB_DOMAIN gismo
      "${_nb_file}"
    )
    target_link_libraries(pygismo_${_src_name} PRIVATE ${PROJECT_NAME} ${_py_target})
    target_link_options(pygismo_${_src_name} PRIVATE ${_py_linker_fix})

    set_target_properties(pygismo_${_src_name} PROPERTIES
      OUTPUT_NAME "${_mod_name}"
      LIBRARY_OUTPUT_DIRECTORY "${PYGISMO_PKG_DIR}"
    )
    file(APPEND "${PYGISMO_PKG_DIR}/__init__.py"
      "from .${_mod_name} import *\n")

    list(APPEND PYGISMO_TARGETS pygismo_${_src_name})
  endforeach()

  file(APPEND "${PYGISMO_PKG_DIR}/__init__.py"
    "import importlib as _importlib, warnings as _warnings\n"
    "class _DeprecatedAlias:\n"
    "    _MAP = {'modelling': 'modeling'}\n"
    "    def __init__(self, mod): self._mod = mod\n"
    "    def __getattr__(self, name):\n"
    "        return getattr(_importlib.import_module('.' + self._mod, __package__), name)\n"
    "def __getattr__(name):\n"
    "    if name in _DeprecatedAlias._MAP:\n"
    "        _warnings.warn(f'pygismo.{name} is deprecated, use pygismo.{_DeprecatedAlias._MAP[name]}', DeprecationWarning, stacklevel=2)\n"
    "        return _DeprecatedAlias(_DeprecatedAlias._MAP[name])\n"
    "    raise AttributeError(f'module pygismo has no attribute {name!r}')\n"
  )

  set(PYGISMO_TARGETS ${PYGISMO_TARGETS} CACHE INTERNAL "nanobind module targets")

  add_custom_target(pygismo_full DEPENDS ${PYGISMO_TARGETS})

endif(GISMO_WITH_NANOBIND)

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

if (GISMO_WITH_ADIFF)
   set(${PROJECT_NAME}_SOURCES ${${PROJECT_NAME}_SOURCES}
     "${gismo_SOURCE_DIR}/external/gsAutoDiff.h")
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
  LABELS "${PROJECT_NAME}"
  FOLDER "G+Smo libraries"
  )
  #generate_export_header(${PROJECT_NAME})



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
      # Note: Download and install "Intel oneAPI Base Toolkit"
      # Then source /path-to/intel/oneapi/setvars.sh
      # use Intel compilers if desired: export CC=icx CXX=icpx
      # Then run e.g.: make ../  -DGISMO_WITH_PARDISO=ON  -DPARDISO_USE_MKL=ON -DEIGEN_USE_MKL_ALL=ON -DCMAKE_BUILD_TYPE=Release -DGISMO_WITH_OPENMP=ON -DMKL_INTERFACE=lp64
      find_package(MKL REQUIRED)
      target_link_libraries(${PROJECT_NAME} MKL::MKL)
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
  # Note: Download and install "Intel oneAPI Base Toolkit"
  # Then source /path-to/intel/oneapi/setvars.sh
  # use Intel compilers if desired: export CC=icx CXX=icpx
  # Then run e.g.: cmake ../  -DGISMO_WITH_PARDISO=ON  -DPARDISO_USE_MKL=ON -DEIGEN_USE_MKL_ALL=ON -DCMAKE_BUILD_TYPE=Release -DGISMO_WITH_OPENMP=ON -DMKL_INTERFACE=lp64
  find_package(MKL REQUIRED)
  target_link_libraries(${PROJECT_NAME} MKL::MKL)
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
