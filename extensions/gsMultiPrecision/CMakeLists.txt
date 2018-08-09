### CMakeLists.txt ---
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2016 - RICAM-Linz.
######################################################################

## Multi-precision floating-point extension
project(gsMultiPrecisionExtension)

# Use internally built GMP/MPFR libraries
set(GISMO_BUILD_MULTIPRECISION BOOL ON CACHE INTERNAL "")

# Collect file names
aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HEADERS)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_SOURCES)
aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HPPFILES)

# Apply same configuration as G+Smo
include(gsConfig)

message("STATUS CMake Will attempt to download GMP sources.")
include(ExternalProject)
ExternalProject_Add(GMP
          SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/gmp
          URL https://gmplib.org/download/gmp/gmp-6.1.2.tar.bz2
          CONFIGURE_COMMAND CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ${CMAKE_CURRENT_BINARY_DIR}/gmp/configure --enable-cxx --prefix=<INSTALL_DIR>
          #BUILD_COMMAND ""
          UPDATE_COMMAND ""
          #INSTALL_COMMAND ""
)

ExternalProject_Get_Property(GMP install_dir)
set(GMP_LIBRARY_DIR ${install_dir}/lib CACHE INTERNAL "")
set(GMP_INCLUDE_DIR ${install_dir}/include CACHE INTERNAL "")
include_directories(${GMP_INCLUDE_DIR})

message("STATUS CMake Will attempt to download MPFR sources.")
include(ExternalProject)
ExternalProject_Add(MPFR
          SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/mpfr
          URL https://www.mpfr.org/mpfr-current/mpfr-4.0.1.tar.gz
          CONFIGURE_COMMAND CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ${CMAKE_CURRENT_BINARY_DIR}/mpfr/configure  --with-gmp-include=${GMP_INCLUDE_DIR} --with-gmp-lib=${GMP_LIBRARY_DIR} --prefix=<INSTALL_DIR>
          #BUILD_COMMAND ""
          UPDATE_COMMAND ""
          #INSTALL_COMMAND ""
)

ExternalProject_Get_Property(MPFR install_dir)
set(MPFR_LIBRARY_DIR ${install_dir}/lib CACHE INTERNAL "")
set(MPFR_INCLUDE_DIR ${install_dir}/include CACHE INTERNAL "")
include_directories(${MPFR_INCLUDE_DIR})

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

set(gismo_EXTENSIONS ${gismo_EXTENSIONS} $<TARGET_OBJECTS:${PROJECT_NAME}>
    CACHE INTERNAL "${PROJECT_NAME} extensions to be included")

  add_dependencies(${PROJECT_NAME} GMP)
  add_dependencies(${PROJECT_NAME} MPFR)
  add_dependencies(MPFR GMP)

install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo/gsMultiPrecision/
        FILES_MATCHING PATTERN "*.h")
