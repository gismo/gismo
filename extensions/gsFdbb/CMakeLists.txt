### CMakeLists.txt ---
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2016 - RICAM-Linz.
######################################################################

## FDBB extension
project(gsFdbbExtension)

# Collect file names
aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HEADERS)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_SOURCES)
aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HPPFILES)

# Apply same configuration as G+Smo
include(gsConfig)

message("STATUS CMake Will attempt to download FDBB sources.")
include(ExternalProject)
ExternalProject_Add(FDBB
          SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/FDBB
          #URL https://gitlab.com/mmoelle1/FDBB/repository/archive.zip
          GIT_REPOSITORY https://gitlab.com/mmoelle1/FDBB.git
          CONFIGURE_COMMAND ""
          BUILD_COMMAND ""
          UPDATE_COMMAND ""
          INSTALL_COMMAND ""
)

set(FDBB_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/FDBB/fdbb CACHE INTERNAL "")
include_directories(${FDBB_INCLUDE_DIR})

add_library(${PROJECT_NAME} OBJECT
    ${${PROJECT_NAME}_SOURCES}
    ${${PROJECT_NAME}_HEADERS}
    ${${PROJECT_NAME}_HPPFILES}
)

set_target_properties(${PROJECT_NAME} PROPERTIES
    COMPILE_DEFINITIONS gismo_EXPORTS 
    POSITION_INDEPENDENT_CODE ON
    LINKER_LANGUAGE CXX
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON
    FOLDER "G+Smo extensions"
    )

if(${CMAKE_CXX_STANDARD} EQUAL "11")
  add_definitions(-DSUPPORT_CXX11)
elseif(${CMAKE_CXX_STANDARD} EQUAL "14")
  add_definitions(-DSUPPORT_CXX14)
else()
  message(FATAL_ERROR "FDDB requires C++11 or C++14 standard")
endif()

set(gismo_EXTENSIONS ${gismo_EXTENSIONS} $<TARGET_OBJECTS:${PROJECT_NAME}>
    CACHE INTERNAL "${PROJECT_NAME} extensions to be included")

add_dependencies(${PROJECT_NAME} FDBB)

install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo/gsFDBB/
        FILES_MATCHING PATTERN "*.h")
