### CMakeLists.txt ---
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2016 - RICAM-Linz.
######################################################################

## Unum extension
project(gsUnumExtension)

# Collect file names
aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HEADERS)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_SOURCES)
aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HPPFILES)

# Apply same configuration as G+Smo
include(gsConfig)

message("STATUS CMake Will attempt to download Unum sources.")
include(ExternalProject)
ExternalProject_Add(Unum
          SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/Unum
          GIT_REPOSITORY https://github.com/stillwater-sc/universal.git
          CONFIGURE_COMMAND ""
          BUILD_COMMAND ""
          UPDATE_COMMAND ""
          INSTALL_COMMAND ""
)

set(UNUM_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/Unum/posit CACHE INTERNAL "")
include_directories(${UNUM_INCLUDE_DIR})

add_library(${PROJECT_NAME} OBJECT
    ${${PROJECT_NAME}_SOURCES}
    ${${PROJECT_NAME}_HEADERS}
    ${${PROJECT_NAME}_HPPFILES}
)

add_dependencies(${PROJECT_NAME} Unum)

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

add_dependencies(${PROJECT_NAME} Unum)

install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo/gsUnum/
        FILES_MATCHING PATTERN "*.h")
