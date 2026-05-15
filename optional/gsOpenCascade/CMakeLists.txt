######################################################################
## CMakeLists.txt --- gsOpenCascade
## This file is part of the G+Smo library.
## 
## Author: Angelos Mantzaflaris 
######################################################################

# Open CASCADE Technology extension
project(gsOpenCascade)

if(CMAKE_CXX_STANDARD EQUAL "98")
  message(FATAL_ERROR "CMAKE_CXX_STANDARD must be 11 or better for OpenCascade extensions to compile.")
endif()

# Apply same configuration as G+Smo
include(gsConfig)

## OpenCascade library
## First try to find a system OpenCASCADE installation.
## We must avoid picking up any OpenCASCADEConfig.cmake that OCC's own cmake
## may have generated inside the gismo build tree (as part of an ExternalProject).
## Strategy: unset any cached OpenCASCADE_DIR, then search ONLY the known system
## locations with NO_DEFAULT_PATH so cmake ignores every stale cache entry and
## the cmake package registry entirely.
set(_occ_system_hints
    /usr/lib/x86_64-linux-gnu/cmake/opencascade
    /usr/lib/cmake/opencascade
    /usr/local/lib/cmake/opencascade
    /usr/share/opencascade/cmake)
unset(OpenCASCADE_DIR CACHE)
find_package(OpenCASCADE QUIET
    NO_DEFAULT_PATH PATHS ${_occ_system_hints})
if (NOT OpenCASCADE_FOUND)
  find_package(OpenCascade QUIET
      NO_DEFAULT_PATH PATHS ${_occ_system_hints})
  if (OpenCascade_FOUND)
    set(OpenCASCADE_FOUND               TRUE)
    set(OpenCASCADE_INCLUDE_DIR         ${OpenCascade_INCLUDE_DIR})
    set(OpenCASCADE_LIBRARY_DIR         ${OpenCascade_LIBRARY_DIR})
    set(OpenCASCADE_MAJOR_VERSION       ${OpenCascade_MAJOR_VERSION})
    set(OpenCASCADE_MINOR_VERSION       ${OpenCascade_MINOR_VERSION})
    set(OpenCASCADE_MAINTENANCE_VERSION ${OpenCascade_MAINTENANCE_VERSION})
  endif()
endif()

if (NOT OpenCASCADE_FOUND)
  # When system gmsh (gsGmsh) is also active, it was compiled against the system OCC.
  # Building OCC 7.3.0 from source will cause an ABI mismatch with newer gmsh versions.
  # Require system OCC dev packages in that case.
  if (DEFINED GMSH_LIB)
    message(FATAL_ERROR
      "System gmsh (${GMSH_LIB}) requires a matching system OpenCascade installation,\n"
      "but no OpenCascade development headers were found on this system.\n"
      "Building OCC 7.3.0 from source would cause an ABI mismatch.\n"
      "On Debian/Ubuntu, install the system OCC dev packages with:\n"
      "  sudo apt-get install libocct-foundation-dev libocct-modeling-algorithms-dev\n"
      "             libocct-modeling-data-dev libocct-data-exchange-dev\n"
      "             libocct-visualization-dev libocct-ocaf-dev\n"
      "Then re-run cmake.")
  endif()

  message(STATUS "System OpenCascade not found. Will download and build from source.")
  # Set OpenCascade version
  set(OCC_VER "7.3.0")

  # Option to control X11/Visualization support
  option(OCC_ENABLE_X11 "Enable X11 support for the OpenCascade Visualization module (requires libx11-dev, libxext-dev, libxmu-dev, libxi-dev)" ON)

  if (OCC_ENABLE_X11)
    find_package(X11 COMPONENTS X11 Xext Xmu Xi QUIET)
    find_package(OpenGL QUIET)
    set(_OCC_MISSING_DEPS "")
    if (NOT X11_FOUND)
      set(_OCC_MISSING_DEPS "${_OCC_MISSING_DEPS} X11 (libx11-dev libxext-dev libxmu-dev libxi-dev)")
    endif()
    if (NOT OpenGL_FOUND OR NOT OpenGL_GLX_FOUND)
      set(_OCC_MISSING_DEPS "${_OCC_MISSING_DEPS} OpenGL/GLX (libgl-dev mesa-common-dev)")
    endif()
    if (_OCC_MISSING_DEPS)
      message(FATAL_ERROR
        "OCC_ENABLE_X11 is ON but the following required development headers are missing:${_OCC_MISSING_DEPS}\n"
        "On Debian/Ubuntu, install them with:\n"
        "  sudo apt-get install libx11-dev libxext-dev libxmu-dev libxi-dev libgl-dev mesa-common-dev\n"
        "On CentOS/Red Hat:\n"
        "  sudo yum install libX11-devel libXext-devel libXmu-devel libXi-devel mesa-libGL-devel\n"
        "Alternatively, set -DOCC_ENABLE_X11=OFF to build without visualization support.")
    endif()
  endif()

  # Download OpenCascade sources at configure time
  include(gsFetch)
  gismo_fetch_directory(OCC
      URL_MD5 "35f3da72c1fb53e32e5cfd826113b32a"
      URL https://files.inria.fr/gismo/gz/opencascade-${OCC_VER}.tgz
      DESTINATION    external
      )
    
    # Set OpenCascade include directory
    set(OCC_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/OCC-prefix/include/opencascade CACHE INTERNAL "")

    # Set OpenCascade libraries
    set(OCC_LIBNAMES TKBO TKBool TKBRep TKernel TKFeat TKFillet TKG2d TKG3d TKGeomAlgo TKGeomBase TKHLR TKIGES TKMath TKMesh TKOffset TKPrim TKShHealing TKSTEP TKSTEPAttr TKSTEPBase TKSTEP209 TKSTL TKTopAlgo TKXSBase TKVRML)
    if (OCC_ENABLE_X11)
      list(APPEND OCC_LIBNAMES TKV3d TKService)
    endif()
    foreach(OCC_LIB ${OCC_LIBNAMES})
      list(APPEND OCC_LIBRARIES ${CMAKE_CURRENT_BINARY_DIR}/OCC-prefix/lib/${CMAKE_SHARED_LIBRARY_PREFIX}${OCC_LIB}${CMAKE_SHARED_LIBRARY_SUFFIX})
    endforeach()

  # Build OpenCascade library at compile time
  include(ExternalProject)
  ExternalProject_Add(OCC
      BINARY_DIR        ${CMAKE_CURRENT_BINARY_DIR}/OCC
      SOURCE_DIR        ${gismo_externals}/OCC
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
                 -DCMAKE_INSTALL_MESSAGE=NEVER
                 -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                 -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                 -DBUILD_MODULE_Draw:BOOL=FALSE
                 -DBUILD_MODULE_Visualization:BOOL=${OCC_ENABLE_X11}
      BUILD_BYPRODUCTS "${OCC_LIBRARIES}"
      DOWNLOAD_COMMAND  ""
      UPDATE_COMMAND    ""
      )

  # Install OpenCascade header files
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/OCC-prefix/include/opencascade
      DESTINATION include/gismo/opencascade/
      FILES_MATCHING PATTERN "*.h")
else()
  message(STATUS "Using system OpenCascade ${OpenCASCADE_MAJOR_VERSION}.${OpenCASCADE_MINOR_VERSION}.${OpenCASCADE_MAINTENANCE_VERSION} from: ${OpenCASCADE_INCLUDE_DIR}")

  # Populate OCC_INCLUDE_DIR and OCC_LIBRARIES from the system installation.
  set(OCC_INCLUDE_DIR ${OpenCASCADE_INCLUDE_DIR} CACHE INTERNAL "")

  # Mirror the same subset of libraries used in the build-from-source path.
  set(OCC_LIBNAMES TKBO TKBool TKBRep TKernel TKFeat TKFillet TKG2d TKG3d TKGeomAlgo TKGeomBase TKHLR TKIGES TKMath TKMesh TKOffset TKPrim TKShHealing TKSTEP TKSTEPAttr TKSTEPBase TKSTEP209 TKSTL TKTopAlgo TKXSBase TKVRML TKV3d TKService)
  foreach(OCC_LIB ${OCC_LIBNAMES})
    list(APPEND OCC_LIBRARIES ${OpenCASCADE_LIBRARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}${OCC_LIB}${CMAKE_SHARED_LIBRARY_SUFFIX})
  endforeach()
  set(OCC_LIBRARIES ${OCC_LIBRARIES} CACHE INTERNAL "")
endif()

##
## gsOpenCascade extension
##
include_directories(SYSTEM ${OCC_INCLUDE_DIR})
set(OCC_SOURCES gsReadOcct.cpp gsWriteOcct.cpp)

# Compile the extension as part of the G+Smo library
add_library(${PROJECT_NAME} OBJECT ${OCC_SOURCES})

if (NOT OpenCASCADE_FOUND)
  add_dependencies(${PROJECT_NAME} OCC)
endif(NOT OpenCASCADE_FOUND)

# Set standard properties for all G+Smo extensions
set_target_properties(${PROJECT_NAME} PROPERTIES
    COMPILE_DEFINITIONS gismo_EXPORTS 
    POSITION_INDEPENDENT_CODE ON
    LINKER_LANGUAGE CXX
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON
    FOLDER "G+Smo extensions" )

#  add_library(TKBO SHARED IMPORTED)
#  add_dependencies(TKBO OCC)
#  set_target_properties(TKBO
#    PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/OCC-prefix/lib/libTKBO.so)

# Add gsOpenCascade extension to the list of G+Smo extensions
set(gismo_EXTENSIONS ${gismo_EXTENSIONS} $<TARGET_OBJECTS:${PROJECT_NAME}>
    CACHE INTERNAL "Gismo extensions to be included")
    
# Add OpenCascade include directory to G+Smo standard include directories
set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${OCC_INCLUDE_DIR}
    CACHE INTERNAL "Gismo include directories")

# Link G+Smo to OpenCascade library (either dynamically or statically)
set(gismo_LINKER ${gismo_LINKER} ${OCC_LIBRARIES}
    CACHE INTERNAL "Gismo extra linker objects")

# Install gsOpenCascade extension header files
install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo/gsOpenCascade/
        FILES_MATCHING PATTERN "*.h")

# Add example files
aux_cpp_directory(${CMAKE_CURRENT_SOURCE_DIR}/examples FILES)
foreach(file ${FILES})
    add_gismo_executable(${file})
    get_filename_component(tarname ${file} NAME_WE)
    set_property(TEST ${tarname} PROPERTY LABELS "${PROJECT_NAME}")
    set_target_properties(${tarname} PROPERTIES FOLDER "${PROJECT_NAME}")
    target_include_directories(${tarname} PRIVATE ${OCC_INCLUDE_DIR})
    target_link_libraries(${tarname} gismo ${OCC_LIBRARIES})
    if (NOT OpenCASCADE_FOUND)
        add_dependencies(${tarname} OCC)
    endif()
    install(TARGETS ${tarname} DESTINATION "${BIN_INSTALL_DIR}" COMPONENT exe OPTIONAL)
endforeach(file ${FILES})
