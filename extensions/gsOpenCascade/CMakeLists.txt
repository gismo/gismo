######################################################################
## CMakeLists.txt --- gsOpenCascade
## This file is part of the G+Smo library.
## 
## Author: Angelos Mantzaflaris 
######################################################################

project(gsOpenCascade)

#find_package(OpenCascade QUIET)

if (NOT OpenCascade_FOUND)
  message("Will attempt to download and compile OpenCascade.")
  include(ExternalProject)
  ExternalProject_Add(OCC
    SOURCE_DIR        ${CMAKE_CURRENT_BINARY_DIR}/OpenCascade
    URL_MD5 "35f3da72c1fb53e32e5cfd826113b32a"
    URL https://people.ricam.oeaw.ac.at/gismo/thirdparty/opencascade-7.3.0.tgz
    UPDATE_COMMAND ""
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> -DCMAKE_INSTALL_MESSAGE=NEVER
               -DOCE_TESTING=OFF -DOCE_VISUALISATION=OFF -DOCE_DISABLE_X11=ON   )

  set(OCC_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/OCC-prefix/include/opencascade CACHE INTERNAL "")
  
  set(OCC_LIBNAMES TKBO TKBool TKBRep TKernel TKFeat TKFillet TKG2d TKG3d TKGeomAlgo TKGeomBase TKHLR TKIGES TKMath TKMesh TKOffset TKPrim TKShHealing TKSTEP TKSTEPAttr TKSTEPBase TKSTEP209 TKSTL TKTopAlgo TKXSBase TKVRML TKV3d TKService )
  foreach(OCC_LIB ${OCC_LIBNAMES})
    list(APPEND OCC_LIBRARIES ${CMAKE_CURRENT_BINARY_DIR}/OCC-prefix/lib/${CMAKE_SHARED_LIBRARY_PREFIX}${OCC_LIB}${CMAKE_SHARED_LIBRARY_SUFFIX})
  endforeach()
  
else()
    #add_custom_target (OCC)
endif()

message("{OCC_INCLUDE_DIR}  ${OCC_INCLUDE_DIR}")
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
    CACHE INTERNAL "Gismo extensions to be included")
set(gismo_LINKER ${gismo_LINKER} ${OCC_LIBRARIES}
    CACHE INTERNAL "${PROJECT_NAME} extra linker objects")

install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo/gsOpenCascade/
        FILES_MATCHING PATTERN "*.h")
