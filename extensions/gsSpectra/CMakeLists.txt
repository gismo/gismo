### CMakeLists.txt --- Spectra library
## Author: Angelos Mantzaflaris 
#
# The library is header-only, the extra code is downloaded in place of
# the source folder of G+Smo
#
######################################################################

## Spectra extension
project(gsSpectraExtension)

# Collect file names
aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HEADERS)
#aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_SOURCES)
#aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HPPFILES)

# Apply same configuration as G+Smo
include(gsConfig)

message(STATUS "CMake Will attempt to download Spectra sources.")
#include(ExternalProject)
gismo_fetch_directory(Spectra
          #SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/Spectra
	  SVN_REPOSITORY https://github.com/yixuan/spectra/trunk/include
	  SVN_TRUST_CERT 1
	  DESTINATION extensions/gsSpectra
)
#endif(NOT SPECTRA_FOUND)


set(SPECTRA_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Spectra/Spectra CACHE INTERNAL "")
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Spectra/Spectra)

add_library(${PROJECT_NAME} INTERFACE)
target_sources(${PROJECT_NAME} INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_HEADERS}>)

target_include_directories(${PROJECT_NAME} INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/Spectra
  #$<INSTALL_INTERFACE:gismo/gsSpectra>
)

add_dependencies(${PROJECT_NAME} Spectra)

#set_target_properties(${PROJECT_NAME} PROPERTIES
#    COMPILE_DEFINITIONS gismo_EXPORTS 
#    POSITION_INDEPENDENT_CODE ON
#    LINKER_LANGUAGE CXX
#    CXX_VISIBILITY_PRESET hidden
#    VISIBILITY_INLINES_HIDDEN ON
#    FOLDER "G+Smo extensions"
#    )

#set(gismo_EXTENSIONS ${gismo_EXTENSIONS} $<TARGET_OBJECTS:${PROJECT_NAME}>
#    CACHE INTERNAL "${PROJECT_NAME} extensions to be included")

set(gismo_LINKER ${gismo_LINKER} ${PROJECT_NAME}
    CACHE INTERNAL "${PROJECT_NAME} extra linker objects")

install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo/gsSpectra
        FILES_MATCHING PATTERN "*.h")
	
