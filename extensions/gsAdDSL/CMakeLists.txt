### CMakeLists.txt ---
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2016 - RICAM-Linz.
######################################################################

if(NOT GISMO_WITH_CODIPACK)
  message(FATAL_ERROR "gsAdDSLExtension cannot be compiled without GISMO_WITH_CODIPACK=ON")
endif()

## AdDSL extension
project(gsAdDSLExtension)

# Collect file names
aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HEADERS)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_SOURCES)
aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HPPFILES)

# Apply same configuration as G+Smo
include(gsConfig)

message("STATUS CMake Will attempt to download AdDSL sources.")
include(ExternalProject)
ExternalProject_Add(AdDSL
          SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/AdDSL
          GIT_REPOSITORY git@max-sagebaum.de:adDsl
          CONFIGURE_COMMAND ""
          BUILD_COMMAND make
          BUILD_IN_SOURCE 1
          UPDATE_COMMAND ""
          INSTALL_COMMAND ""
)

set(ADDSL_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/AdDSL/include CACHE INTERNAL "")
include_directories(${ADDSL_INCLUDE_DIR})

add_library(${PROJECT_NAME} OBJECT
    ${${PROJECT_NAME}_SOURCES}
    ${${PROJECT_NAME}_HEADERS}
    ${${PROJECT_NAME}_HPPFILES}
)

add_dependencies(${PROJECT_NAME} AdDSL)
add_dependencies(AdDSL gsCoDiPackExtension)

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
        DESTINATION include/gismo/gsAdDSL/
        FILES_MATCHING PATTERN "*.h")
