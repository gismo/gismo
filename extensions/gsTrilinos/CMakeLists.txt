### CMakeLists.txt ---
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2016 - RICAM-Linz.
######################################################################

cmake_minimum_required(VERSION 2.8.10.2)

if(POLICY CMP0048)# CMake 3.0
cmake_policy(SET CMP0011 OLD)
cmake_policy(SET CMP0042 OLD)
cmake_policy(SET CMP0048 OLD)
endif()

if(POLICY CMP0054)# CMake 3.1
cmake_policy(SET CMP0054 OLD)
endif()

if(POLICY CMP0053)# CMake 3.1.3
cmake_policy(SET CMP0053 OLD)
endif()

if(POLICY CMP0063)# CMake 3.3
cmake_policy(SET CMP0063 NEW)
endif()

## Trilinos extension
project(gsTrilinosExtension)

# Collect file names
aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HEADERS)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_SOURCES)
aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HPPFILES)

# Apply same configuration as G+Smo
include( gsConfig)

# Look for MPI (now: handled by GISMO_WITH_MPI)
#find_package(MPI REQUIRED)
#include_directories(SYSTEM ${MPI_INCLUDE_PATH})

# Look for trilinos library
#find_package(TRILINOS REQUIRED) # for testing user-provided mode
find_package(TRILINOS QUIET)
#find_package(BLAS REQUIRED)   # needs .._DIR on MS windows
#find_package(LAPACK REQUIRED) # needs LAPACK_DIR on MS windows

#
find_path(TPL_MPI_INCLUDE_DIRS mpi.h ${MPI_INCLUDE_PATH})
list(GET MPI_CXX_LIBRARIES 0 TPL_MPI_LIBRARY_DIRS)
get_filename_component(TPL_MPI_LIBRARY_DIRS ${TPL_MPI_LIBRARY_DIRS} DIRECTORY)
get_filename_component(TPL_MPI_BASE_DIR ${TPL_MPI_LIBRARY_DIRS} DIRECTORY)

if (NOT GISMO_BUILD_CPP11)
message(WARNING "GISMO_BUILD_CPP11 must be ON for Trilinos extensions to compile ?.")
endif()

# Note: might need export FC=/usr/bin/gfortran
#enable_language(Fortran)

if (NOT TRILINOS_FOUND)
set(TRILINOS_VER "12-8-1")
set(TRILINOS_MD5 "01c0026f1e2050842857db941060ecd5")
message(STATUS "CMake Will attempt to download and compile Trilinos ${TRILINOS_VER}.")
include(ExternalProject)
ExternalProject_Add(Trilinos
          SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/Trilinos
          URL_MD5 ${TRILINOS_MD5}
          URL https://github.com/trilinos/Trilinos/archive/trilinos-release-${TRILINOS_VER}.tar.gz
          CMAKE_ARGS   -DCMAKE_BUILD_TYPE=Release   
                       -DTrilinos_ENABLE_TESTS=OFF  
                       -DTrilinos_ENABLE_Epetra=ON
                       -DTrilinos_ENABLE_Amesos=ON
                       -DTrilinos_ENABLE_AztecOO=ON
                       -DTrilinos_ENABLE_Anasazi=ON  
                       -DTrilinos_ENABLE_Belos=${GISMO_BUILD_CPP11}
                       -DTrilinos_ENABLE_NOX=OFF
                       #-DTrilinos_ENABLE_ML=ON # multigrid
                       #-DTPL_ENABLE_Pthread=OFF # Note: MS Win might not have Pthread
                       -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
                       -DBUILD_SHARED_LIBS:BOOL=ON
                       -DTPL_ENABLE_MPI=${GISMO_WITH_MPI}
                       -DTPL_MPI_INCLUDE_DIRS=${TPL_MPI_INCLUDE_DIRS}
                       -DTPL_MPI_LIBRARY_DIRS=${TPL_MPI_LIBRARY_DIRS}
                       #-DMPI_BASE_DIR=${TPL_MPI_BASE_DIR}
                       #-DTPL_MPI_LIBRARIES=${MPI_CXX_LIBRARIES} # bug ? 2 elements
                       -DTrilinos_ENABLE_CXX11=${GISMO_BUILD_CPP11} # Note: this must be ON, otherwise Trilinos compilation fails. Must be OFF on MS Windows
		       -DTPL_ENABLE_SuperLU=${GISMO_WITH_SUPERLU}
	 	       ## Point to specific paths (! no space after "=")
		       -DTPL_LAPACK_LIBRARIES=${LAPACK_LIBRARIES}
      		       -DTPL_BLAS_LIBRARIES=${LAPACK_LIBRARIES}
		       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
		       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
		       #-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
          UPDATE_COMMAND ""
	  #DOWNLOAD_NO_PROGRESS 0
)
      set(TRILINOS_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/Trilinos-prefix/include CACHE INTERNAL "")
      set(TRILINOS_LIBRARIES ${CMAKE_CURRENT_BINARY_DIR}/Trilinos-prefix/lib/${CMAKE_SHARED_LIBRARY_PREFIX}epetra${CMAKE_SHARED_LIBRARY_SUFFIX};
      ${CMAKE_CURRENT_BINARY_DIR}/Trilinos-prefix/lib/${CMAKE_SHARED_LIBRARY_PREFIX}teuchoscore${CMAKE_SHARED_LIBRARY_SUFFIX};
      ${CMAKE_CURRENT_BINARY_DIR}/Trilinos-prefix/lib/${CMAKE_SHARED_LIBRARY_PREFIX}teuchosnumerics${CMAKE_SHARED_LIBRARY_SUFFIX};
      ${CMAKE_CURRENT_BINARY_DIR}/Trilinos-prefix/lib/${CMAKE_SHARED_LIBRARY_PREFIX}kokkoscore${CMAKE_SHARED_LIBRARY_SUFFIX};
      ${CMAKE_CURRENT_BINARY_DIR}/Trilinos-prefix/lib/${CMAKE_SHARED_LIBRARY_PREFIX}teuchoskokkoscomm${CMAKE_SHARED_LIBRARY_SUFFIX};
      ${CMAKE_CURRENT_BINARY_DIR}/Trilinos-prefix/lib/${CMAKE_SHARED_LIBRARY_PREFIX}teuchosparameterlist${CMAKE_SHARED_LIBRARY_SUFFIX};
      ${CMAKE_CURRENT_BINARY_DIR}/Trilinos-prefix/lib/${CMAKE_SHARED_LIBRARY_PREFIX}teuchoscomm${CMAKE_SHARED_LIBRARY_SUFFIX};
      ${CMAKE_CURRENT_BINARY_DIR}/Trilinos-prefix/lib/${CMAKE_SHARED_LIBRARY_PREFIX}teuchoskokkoscompat${CMAKE_SHARED_LIBRARY_SUFFIX};
      ${CMAKE_CURRENT_BINARY_DIR}/Trilinos-prefix/lib/${CMAKE_SHARED_LIBRARY_PREFIX}teuchosremainder${CMAKE_SHARED_LIBRARY_SUFFIX}
      ${CMAKE_CURRENT_BINARY_DIR}/Trilinos-prefix/lib/${CMAKE_SHARED_LIBRARY_PREFIX}amesos${CMAKE_SHARED_LIBRARY_SUFFIX}
      ${CMAKE_CURRENT_BINARY_DIR}/Trilinos-prefix/lib/${CMAKE_SHARED_LIBRARY_PREFIX}aztecoo${CMAKE_SHARED_LIBRARY_SUFFIX}
      ${CMAKE_CURRENT_BINARY_DIR}/Trilinos-prefix/lib/${CMAKE_SHARED_LIBRARY_PREFIX}epetraext${CMAKE_SHARED_LIBRARY_SUFFIX}
      ${CMAKE_CURRENT_BINARY_DIR}/Trilinos-prefix/lib/${CMAKE_SHARED_LIBRARY_PREFIX}belos${CMAKE_SHARED_LIBRARY_SUFFIX}
      ${CMAKE_CURRENT_BINARY_DIR}/Trilinos-prefix/lib/${CMAKE_SHARED_LIBRARY_PREFIX}anasazi${CMAKE_SHARED_LIBRARY_SUFFIX}
      CACHE INTERNAL "")

#      set(TRILINOS_DIR ${CMAKE_CURRENT_BINARY_DIR}/Trilinos-prefix CACHE INTERNAL "")

      #add_library(${LIBNAME} SHARED IMPORTED )
      #set_target_properties(${LIBNAME} PROPERTIES IMPORTED_LOCATION ${install_dir}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}${LIBNAME}${CMAKE_SHARED_LIBRARY_SUFFIX} )
  
#else()
endif(NOT TRILINOS_FOUND)

#include_directories(${EPETRA_INCLUDE_PATH})
include_directories(${TRILINOS_INCLUDE_DIR})

# Compile statically the extension
add_library(${PROJECT_NAME} OBJECT
${${PROJECT_NAME}_SOURCES}
${${PROJECT_NAME}_HEADERS}
${${PROJECT_NAME}_HPPFILES}
#${TRILINOS_LIBRARIES}
)

if (NOT TRILINOS_FOUND)
  add_dependencies(${PROJECT_NAME} Trilinos)
  set(TRILINOS_FOUND TRUE)
endif(NOT TRILINOS_FOUND)

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
set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS}
  CACHE INTERNAL "${PROJECT_NAME} include directories")

set(gismo_LINKER ${gismo_LINKER} ${LAPACK_LIBRARIES} ${TRILINOS_LIBRARIES}
    CACHE INTERNAL "${PROJECT_NAME} extra linker objects")

install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo/gsTrilinos/
        FILES_MATCHING PATTERN "*.h")

#set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

#message("TRILINOS_LIBRARIES ${TRILINOS_LIBRARIES}")
#list(LENGTH TRILINOS_LIBRARIES tlen)
#message("NUM TRILINOS_LIBRARIES ${tlen}")
