######################################################################
### CMakeLists.txt --- gsTrilinos
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
######################################################################

cmake_minimum_required(VERSION 2.8.12)

if(POLICY CMP0048)# CMake 3.0
cmake_policy(SET CMP0011 NEW)
cmake_policy(SET CMP0042 NEW)
cmake_policy(SET CMP0048 NEW)
endif()

if(POLICY CMP0054)# CMake 3.1
cmake_policy(SET CMP0054 NEW)
endif()

if(POLICY CMP0053)# CMake 3.1.3
cmake_policy(SET CMP0053 NEW)
endif()

if(POLICY CMP0063)# CMake 3.3
cmake_policy(SET CMP0063 NEW)
endif()

if(POLICY CMP0077)
cmake_policy(SET CMP0077 NEW)
endif()

## Trilinos extension
project(gsTrilinosExtension)

# Collect file names
aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HEADERS)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_SOURCES)
aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HPPFILES)

# Apply same configuration as G+Smo
include(gsConfig)

if(NOT GISMO_WITH_MPI)
  message(FATAL_ERROR "Trilinos requires GISMO_WITH_MPI=ON. Please enable MPI and try again.")
endif()

# Look for pre-installed Trilinos libraries
find_package(Trilinos QUIET)

##
## Trilinos packages
##
if (NOT Trilinos_FOUND)

  if(CMAKE_CXX_STANDARD EQUAL "98")
    message(FATAL_ERROR "CMAKE_CXX_STANDARD must be 11 or better for Trilinos extensions to compile.")
  endif()

  # Set Trilinos version
  set(Trilinos_VER "13-0-1")
  set(Trilinos_MD5 "9d76b494fe1c00cedb93d51dfc666a27") #"de912cca43c2ca3b74aa08528ac39dbd")

  # Download Trilinos sources at configure time
  include(gsFetch)
  gismo_fetch_directory(Trilinos
    URL            https://github.com/trilinos/Trilinos/archive/trilinos-release-${Trilinos_VER}.tar.gz
    URL_MD5        ${Trilinos_MD5}
    DESTINATION    external
    PATCH_COMMAND  patch -p1 -N -d ${gismo_externals}/Trilinos < ${CMAKE_CURRENT_SOURCE_DIR}/Trilinos.patch
    )

  # Global Trilinos configuration
  set(Trilinos_BUILD_SHARED_LIBS                 OFF)
  set(Trilinos_ENABLE_TESTS                      OFF)

  if(CMAKE_CXX_STANDARD     EQUAL 11)
    set(Trilinos_ENABLE_CXX11                     ON)
  elseif(CMAKE_CXX_STANDARD EQUAL 14)
    set(Trilinos_ENABLE_CXX14                     ON)
  elseif(CMAKE_CXX_STANDARD EQUAL 17)
    set(Trilinos_ENABLE_CXX17                     ON)
  elseif(CMAKE_CXX_STANDARD EQUAL 20)
    set(Trilinos_ENABLE_CXX20                     ON)
  else()
    message(FATAL_ERROR "Unsupported C++ standard.")
  endif()    

  # Smart reference-counted pointer class
  set(Trilinos_ENABLE_Teuchos                     ON)
  #set(Trilinos_ENABLE_TeuchosParser               ON)
  
  # Linear algebra backends
  set(Trilinos_ENABLE_Epetra                      ON) # non-templated version
  set(Trilinos_ENABLE_Tpetra                      ON) # templated version
  #set(Trilinos_ENABLE_COMPLEX_DOUBLE              ON) # needed by Anasazi
  #set(Trilinos_ENABLE_COMPLEX_FLOAT               ON) # needed by Anasazi
  #set(Teuchos_ENABLE_COMPLEX                      ON)
  #set(Teuchos_ENABLE_DOUBLE                       ON)
  #set(Teuchos_ENABLE_FLOAT                        ON)
  #set(Tpetra_INST_COMPLEX_DOUBLE                  ON)
  #set(Tpetra_INST_COMPLEX_FLOAT                   ON)
  #set(Tpetra_INST_DOUBLE                          ON)
  #set(Tpetra_INST_FLOAT                           ON)
  
  # Sparse direct solvers
  set(Trilinos_ENABLE_Amesos                      ON)
  
  # Iterative solvers
  set(Trilinos_ENABLE_AztecOO                     ON) # Krylov subspace solvers
  set(Trilinos_ENABLE_Belos                       ON) # Krylov subspace solvers
  set(Trilinos_ENABLE_ML                          ON) # Multilevel solver
  set(Trilinos_ENABLE_Ifpack                      ON) # Ifpack
  
  # Eigenvalue solvers
  set(Trilinos_ENABLE_Anasazi                     ON)
  
  # Nonlinear solvers
  set(Trilinos_ENABLE_NOX                         ON)
  set(NOX_ENABLE_ThyraCore                        ON)
  set(NOX_ENABLE_ThyraTpetraAdapters              ON)
  set(NOX_ENABLE_ThyraEpetraExtAdapters           ON)
  
  # Disable all optional packages by default
  set(Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES      OFF)

  # OpenMP
  set(TPL_ENABLE_OpenMP                          ${GISMO_WITH_OPENMP})
  
  # MPI
  #find_path(TPL_MPI_INCLUDE_DIRS mpi.h ${MPI_INCLUDE_PATH})
  #list(GET MPI_CXX_LIBRARIES 0 TPL_MPI_LIBRARY_DIRS)
  #get_filename_component(TPL_MPI_LIBRARY_DIRS ${TPL_MPI_LIBRARY_DIRS} PATH)

  set(TPL_ENABLE_MPI                              ${GISMO_WITH_MPI})
  #set(TPL_MPI_INCLUDE_DIRS                        ${TPL_MPI_INCLUDE_DIRS})
  #set(TPL_MPI_LIBRARY_DIRS                        ${TPL_MPI_LIBRARY_DIRS})
  
  ## Point to specific paths (! no space after "=")
  set(TPL_LAPACK_LIBRARIES                        ${LAPACK_LIBRARIES})
  set(TPL_BLAS_LIBRARIES                          ${LAPACK_LIBRARIES})
  
  # Pardiso solvers
  set(TPL_ENABLE_PARDISO                          ${GISMO_WITH_PARDISO})
  
  # SuperLU solvers
  set(TPL_ENABLE_SuperLU                          ${GISMO_WITH_SUPERLU})
  set(TPL_SuperLU_INCLUDE_DIRS                    ${SUPERLU_INCLUDES})
  set(TPL_SuperLU_LIBRARIES                       ${SUPERLU_LIBRARIES})
    
  # Taucs solvers
  set(TPL_ENABLE_TAUCS                            ${GISMO_WITH_TAUCS})

  # Umfpack solvers
  set(TPL_ENABLE_UMFPACK                          ${GISMO_WITH_UMFPACK})

  # Build Trilinos at compile time
  add_subdirectory(${gismo_externals}/Trilinos ${gismo_BINARY_DIR}/external/Trilinos)
  
  foreach(target ${Trilinos_LIBRARIES})
    set_property(TARGET ${target} PROPERTY POSITION_INDEPENDENT_CODE ON)
  endforeach()

endif(NOT Trilinos_FOUND)

include_directories(${Trilinos_INCLUDE_DIRS})

##
## gsTrilinos extension
##

# Compile gsTrilinos extension as part of the G+Smo library
add_library(${PROJECT_NAME} OBJECT
    ${${PROJECT_NAME}_SOURCES}
    ${${PROJECT_NAME}_HEADERS}
    ${${PROJECT_NAME}_HPPFILES}
    )

# Set standard properties for all G+Smo extensions
set_target_properties(${PROJECT_NAME} PROPERTIES
    COMPILE_DEFINITIONS gismo_EXPORTS
    POSITION_INDEPENDENT_CODE ON
    LINKER_LANGUAGE CXX
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON
    FOLDER "G+Smo extensions"
    )

# Add gsTrilinos extension to the list of G+Smo extensions
set(gismo_EXTENSIONS ${gismo_EXTENSIONS} $<TARGET_OBJECTS:${PROJECT_NAME}>
    CACHE INTERNAL "gismo extensions to be included")

# Add Trilinos include directories to G+Smo standard include directories
set(GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${Trilinos_INCLUDE_DIRS}
    CACHE INTERNAL "gismo include directories")

# Link G+Smo to Trilinos libraries (either dynamically or statically)
set(gismo_LINKER ${gismo_LINKER} ${LAPACK_LIBRARIES} ${Trilinos_LIBRARIES}
    CACHE INTERNAL "gismo extra linker objects")

# Install gsTrilinos extension header files
install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo/gsTrilinos/
        FILES_MATCHING PATTERN "*.h")

##
## Add dependencies
##
if(GISMO_WITH_TRILINOS AND NOT Trilinos_FOUND)
  add_dependencies(${PROJECT_NAME} ${Trilinos_LIBRARIES})
endif(GISMO_WITH_TRILINOS AND NOT Trilinos_FOUND)
