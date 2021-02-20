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

# Install CoDiPack header files
install(DIRECTORY ${CODIPACK_INCLUDE_DIR}/
        DESTINATION include/gismo/ )

# Install gsCoDiPack header file
install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo/
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")

# Add filedata folder
add_definitions(-DELAST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/filedata/") 

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
