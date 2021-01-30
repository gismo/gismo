######################################################################
### CMakeLists.txt --- gsCodiPack
## This file is part of the G+Smo library.
## 
## Author: Angelos Mantzaflaris 
######################################################################

## CoDiPack extension
project(gsCoDiPackExtension)

if(CMAKE_CXX_STANDARD EQUAL "98" )
  message(FATAL_ERROR "CMAKE_CXX_STANDARD must be 11 or better for CoDiPack extensions to compile.")
endif()

# Apply same configuration as G+Smo
include(gsConfig)
set(CMAKE_CXX_VISIBILITY_PRESET default)# this is needed to get the CoDiPack tapes exported

##
## CoDiPack package
##

# Download CoDiPack sources at configure time
include(gsFetch)
gismo_fetch_directory(CoDiPack
    URL https://github.com/SciCompKL/CoDiPack/archive/numericExtension.zip
    DESTINATION    external
)

# Set CoDiPack include directory
set(CODIPACK_INCLUDE_DIR ${gismo_SOURCE_DIR}/external/CoDiPack/include CACHE INTERNAL "")

# Add CoDiPack include directory to G+Smo standard include directories
set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${CODIPACK_INCLUDE_DIR}
  CACHE INTERNAL "Gismo include directories")

# Add needed instantiations to the precompiled library
set (GISMO_EXTRA_INSTANCE ${GISMO_EXTRA_INSTANCE} "codi::RealForward;codi::RealReverse"
  CACHE INTERNAL "Additional instantiations")

add_library(${PROJECT_NAME} OBJECT gsCoDiPack_.cpp)

set_target_properties(${PROJECT_NAME} PROPERTIES
  COMPILE_DEFINITIONS gismo_EXPORTS
  POSITION_INDEPENDENT_CODE ON
  LINKER_LANGUAGE CXX
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN ON
  FOLDER "G+Smo extensions"  )

set(gismo_EXTENSIONS ${gismo_EXTENSIONS} $<TARGET_OBJECTS:${PROJECT_NAME}>
  CACHE INTERNAL "Gismo extensions to be included")

# Install CoDiPack header files
install(DIRECTORY ${CODIPACK_INCLUDE_DIR}/
        DESTINATION include/gismo/ )

# Install gsCoDiPack header file
install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo/
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")

# Add example files
include_directories(${CODIPACK_INCLUDE_DIR})
aux_cpp_directory(${CMAKE_CURRENT_SOURCE_DIR}/examples FILES)
foreach(file ${FILES})
  add_gismo_executable(${file})
  get_filename_component(tarname ${file} NAME_WE) # name without extension
  set_property(TEST ${tarname} PROPERTY LABELS "${PROJECT_NAME}")
  set_target_properties(${tarname} PROPERTIES FOLDER "${PROJECT_NAME}")
  # Install the example executables (optionally)
  install(TARGETS ${tarname} DESTINATION "${BIN_INSTALL_DIR}" COMPONENT exe OPTIONAL)
endforeach(file ${FILES})
