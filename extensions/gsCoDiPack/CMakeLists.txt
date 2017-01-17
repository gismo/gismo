### CMakeLists.txt ---
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2016 - RICAM-Linz.
######################################################################

## CoDiPack extension
project(gsCoDiPackExtension)

# Collect file names
aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HEADERS)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_SOURCES)
aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HPPFILES)

# Apply same configuration as G+Smo
include(gsConfig)

# Look for CoDiPack library
#find_package(CODIPACK QUIET)
#find_package(LAPACK REQUIRED)

#if (NOT CODIPACK_FOUND)

message("STATUS CMake Will attempt to download CoDiPack sources.")
include(ExternalProject)
ExternalProject_Add(CoDiPack
          SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/CoDiPack
          URL https://github.com/SciCompKL/CoDiPack/archive/numericExtension.zip
          CONFIGURE_COMMAND ""
          BUILD_COMMAND ""
          UPDATE_COMMAND ""
          INSTALL_COMMAND ""
)
#endif(NOT CODIPACK_FOUND)

set(CODIPACK_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/CoDiPack/include CACHE INTERNAL "")
include_directories(${CODIPACK_INCLUDE_DIR})

add_library(${PROJECT_NAME} OBJECT
    ${${PROJECT_NAME}_SOURCES}
    ${${PROJECT_NAME}_HEADERS}
    ${${PROJECT_NAME}_HPPFILES}
)

add_dependencies(${PROJECT_NAME} CoDiPack)

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
        DESTINATION include/gismo/gsCoDiPack/
        FILES_MATCHING PATTERN "*.h")
