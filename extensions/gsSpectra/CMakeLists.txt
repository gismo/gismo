### CMakeLists.txt ---
## Author: Angelos Mantzaflaris 
## Copyright (C) 2018 - RICAM-Linz.
######################################################################

## Spectra extension
project(gsSpectraExtension)

# Collect file names
aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HEADERS)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_SOURCES)
#aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HPPFILES)

# Apply same configuration as G+Smo
include(gsConfig)

message("STATUS CMake Will attempt to download Spectra sources.")
include(ExternalProject)
ExternalProject_Add(Spectra
          SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/Spectra
	  SVN_REPOSITORY https://github.com/yixuan/spectra/trunk/include
          CONFIGURE_COMMAND ""
          BUILD_COMMAND ""
          UPDATE_COMMAND ""
          INSTALL_COMMAND ""
)
#endif(NOT SPECTRA_FOUND)

set(SPECTRA_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/Spectra CACHE INTERNAL "")
include_directories(${SPECTRA_INCLUDE_DIR})

add_library(${PROJECT_NAME} OBJECT
    ${${PROJECT_NAME}_HEADERS}
    ${${PROJECT_NAME}_SOURCES} #need source to enforce dependancy
#    ${${PROJECT_NAME}_HPPFILES}
)

add_dependencies(${PROJECT_NAME} Spectra)

set_target_properties(${PROJECT_NAME} PROPERTIES
    COMPILE_DEFINITIONS gismo_EXPORTS 
    POSITION_INDEPENDENT_CODE ON
    LINKER_LANGUAGE CXX
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON
    FOLDER "G+Smo extensions"
    )

set(gismo_EXTENSIONS ${gismo_EXTENSIONS} $<TARGET_OBJECTS:${PROJECT_NAME}>
    CACHE INTERNAL "${PROJECT_NAME} extensions to be included")

install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo/gsSpectra/
        FILES_MATCHING PATTERN "*.h")
