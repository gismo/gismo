######################################################################
## CMakeLists.txt --- gsOpenCascade
## This file is part of the G+Smo library.
## 
## Author: Angelos Mantzaflaris 
######################################################################

project(gsOpenCascade)

if(CMAKE_CXX_STANDARD EQUAL "98")
  message(FATAL_ERROR "CMAKE_CXX_STANDARD must be 11 or better for OpenCascade extensions to compile.")
endif()

# Apply same configuration as G+Smo
include(gsConfig)

##
## OpenCascade library
##
find_package(OpenCascade 7.3.0 QUIET)

if (NOT OpenCascade_FOUND)
  # Set OpenCascade version
  set(OCC_VER "7.3.0")

  # Download OpenCascade sources at configure time
  include(gsFetch)
  gismo_fetch_directory(OCC
      URL_MD5 "35f3da72c1fb53e32e5cfd826113b32a"
      URL https://people.ricam.oeaw.ac.at/gismo/thirdparty/opencascade-${OCC_VER}.tgz
      DESTINATION    external
      )

  # Build OpenCascade library at compile time
  include(ExternalProject)
  ExternalProject_Add(OCC
      BINARY_DIR        ${CMAKE_CURRENT_BINARY_DIR}/OCC
      SOURCE_DIR        ${gismo_externals}/OCC    
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
                 -DCMAKE_INSTALL_MESSAGE=NEVER
                 -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                 -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      DOWNLOAD_COMMAND     ""
      UPDATE_COMMAND ""
      )

  # Set OpenCascade include directory
  set(OCC_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/OCC-prefix/src/OCC-build/include/opencascade CACHE INTERNAL "")

  # Set OpenCascade libraries
  set(OCC_LIBNAMES TKBO TKBool TKBRep TKernel TKFeat TKFillet TKG2d TKG3d TKGeomAlgo TKGeomBase TKHLR TKIGES TKMath TKMesh TKOffset TKPrim TKShHealing TKSTEP TKSTEPAttr TKSTEPBase TKSTEP209 TKSTL TKTopAlgo TKXSBase TKVRML TKV3d TKService )
  foreach(OCC_LIB ${OCC_LIBNAMES})
    list(APPEND OCC_LIBRARIES ${CMAKE_CURRENT_BINARY_DIR}/OCC-prefix/lib/${CMAKE_SHARED_LIBRARY_PREFIX}${OCC_LIB}${CMAKE_SHARED_LIBRARY_SUFFIX})
  endforeach()
endif()

include_directories(SYSTEM ${OCC_INCLUDE_DIR})

set(OCC_SOURCES gsReadBrep.cpp)

add_library(${PROJECT_NAME} OBJECT ${OCC_SOURCES})

if (NOT OpenCascade_FOUND)
  add_dependencies(${PROJECT_NAME} OCC)
endif(NOT OpenCascade_FOUND)

set_target_properties(${PROJECT_NAME} PROPERTIES
    COMPILE_DEFINITIONS gismo_EXPORTS 
    POSITION_INDEPENDENT_CODE ON
    LINKER_LANGUAGE CXX
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON
    FOLDER "G+Smo extensions"  )

# Apply same configuration as G+Smo
include( gsConfig)

set(gismo_EXTENSIONS ${gismo_EXTENSIONS} $<TARGET_OBJECTS:${PROJECT_NAME}>
    CACHE INTERNAL "gismo extensions to be included")
set(gismo_LINKER ${gismo_LINKER} ${OCC_LIBRARIES}
    CACHE INTERNAL "gismo extra linker objects")

install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo/gsOpenCascade/
        FILES_MATCHING PATTERN "*.h")
