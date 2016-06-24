### CMakeLists.txt ---
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2016 - RICAM-Linz.
######################################################################

## Trilinos extension
project(gsTrilinosExtension)

# Collect file names
aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HEADERS)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_SOURCES)
aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HPPFILES)


# Look for trilinos library
find_package(TRILINOS QUIET)
find_package(LAPACK REQUIRED)

if (NOT TRILINOS_FOUND)
  set(TRILINOS_VER "12-6-3")
  message("Trilinos not found (define TRILINOS_DIR if it exists on the system).")
  message("Will attempt to download and compile Trilinos ${TRILINOS_VER} (do: make Trilinos; cmake .")
include(ExternalProject)
ExternalProject_Add(Trilinos
          SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/Trilinos
          URL https://github.com/trilinos/Trilinos/archive/trilinos-release-${TRILINOS_VER}.tar.gz
	  CMAKE_ARGS   -DCMAKE_BUILD_TYPE=Release   
	               -DBUILD_SHARED_LIBS:BOOL=ON  
                       -DTrilinos_ENABLE_TESTS=OFF  
                       -DTrilinos_ENABLE_Epetra=ON
                       #-DCMAKE_INSTALL_PREFIX=
                       #-DTrilinos_ENABLE_CXX11=
          UPDATE_COMMAND ""
)
      set(TRILINOS_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/Trilinos-prefix/include CACHE INTERNAL "")
      set(TRILINOS_LIBRARIES ${CMAKE_BINARY_DIR}/lib/libtrilinos.so  CACHE INTERNAL "")
      
      add_dependencies(${PROJECT_NAME} Trilinos)
#else()
#
endif(NOT TRILINOS_FOUND)

#include_directories(${EPETRA_INCLUDE_PATH})
include_directories(${TRILINOS_INCLUDE_DIR})

# Apply same configuration as G+Smo
include( gsConfig)

# Compile statically the extension
add_library(${PROJECT_NAME} OBJECT
${${PROJECT_NAME}_SOURCES}
${${PROJECT_NAME}_HEADERS}
${${PROJECT_NAME}_HPPFILES}
#${TRILINOS_LIBRARIES}
)

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
set(gismo_LINKER ${gismo_LINKER} ${TRILINOS_LIBRARIES} ${LAPACK_LIBRARIES}
    CACHE INTERNAL "${PROJECT_NAME} extra linker objects")

install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo/gsTrilinos/
        FILES_MATCHING PATTERN "*.h")

#set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

#message("TRILINOS_LIBRARIES ${TRILINOS_LIBRARIES}")
#list(LENGTH TRILINOS_LIBRARIES tlen)
#message("NUM TRILINOS_LIBRARIES ${tlen}")
