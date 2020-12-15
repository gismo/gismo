######################################################################
### CMakeLists.txt --- gsCoreLibs
## This file is part of the G+Smo library.
## 
## Author: Angelos Mantzaflaris, Matthias Moller
######################################################################

# Core libraries extension
project(gsCoreLibsExtension)

# Collect file names
aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HEADERS)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_SOURCES)
aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HPPFILES)

# Apply same configuration as G+Smo
include(gsConfig)

##
## BLAS/Lapack library
##
#if(GISMO_WITH_SUPERLU OR GISMO_WITH_TAUCS OR GISMO_WITH_UMFPACK)
find_package(BLAS QUIET)
find_package(Lapack QUIET)


  if(BLAS_FOUND OR Lapack_FOUND)
    # Set OpenBLAS version
    set(OPENBLAS_VER "0.3.13")

    # Download OpenBLAS sources at configure time
    include(gsFetch)
    gismo_fetch_directory(openblas
        URL            https://github.com/xianyi/OpenBLAS/releases/download/v${OPENBLAS_VER}/OpenBLAS-${OPENBLAS_VER}.tar.gz
        DESTINATION    external
        )

    # Build OpenBLAS library at compile time
    add_subdirectory(${gismo_externals}/openblas ${gismo_BINARY_DIR}/external/openblas)
      
  endif(BLAS_FOUND OR Lapack_FOUND)

#endif()
  
##
## Umfpack library
##
if(GISMO_WITH_UMFPACK)
  # Look for pre-installed Umfpack library
  find_package(Umfpack QUIET)
  
  if(NOT Umfpack_FOUND)
    # Set Umfpack version
    set(UMFPACK_VER "4.5.6")  
    
    # Download Umfpack sources at configure time
    include(gsFetch)
    gismo_fetch_directory(umfpack
        URL            https://people.engr.tamu.edu/davis/SuiteSparse/SuiteSparse-${UMFPACK_VER}.tar.gz
        DESTINATION    external
        )

    # Build Umfpack library at compile time
    include(ExternalProject)
    ExternalProject_Add(umfpack
        BINARY_DIR           ${CMAKE_CURRENT_BINARY_DIR}/umfpack
        SOURCE_DIR           ${gismo_externals}/umfpack
        CONFIGURE_COMMAND    CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ${gismo_externals}/gmp/configure --enable-cxx --enable-shared=no --enable-static=yes --with-pic --prefix=<INSTALL_DIR>
        DOWNLOAD_COMMAND     ""
        UPDATE_COMMAND       ""
        )
    
    # Set Umfpack libraries
    set(UMFPACK_LIBRARY ${CMAKE_CURRENT_BINARY_DIR}/umfpack-prefix/lib/${CMAKE_STATIC_LIBRARY_PREFIX}umfpack${CMAKE_STATIC_LIBRARY_SUFFIX} CACHE INTERNAL "")
    
    # Set Umfpack library and include directories
    ExternalProject_Get_Property(umfpack install_dir)
    set(UMFPACK_LIBRARY_DIR ${install_dir}/lib CACHE INTERNAL "")
    set(UMFPACK_INCLUDE_DIR ${install_dir}/include CACHE INTERNAL "")
    include_directories(${UMFPACK_INCLUDE_DIR})
    
    # Install UMFPACK header files
    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/umfpack-prefix/include
        DESTINATION include/gismo/
        FILES_MATCHING PATTERN "*.h")
  endif(NOT Umfpack_FOUND)
endif(GISMO_WITH_UMFPACK)

##
## gsCoreLibs extension
##

# Compile gsCoreLibs extension as part of the G+Smo library
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

# Add gsCoreLibs extension to the list of G+Smo extensions
set(gismo_EXTENSIONS ${gismo_EXTENSIONS} $<TARGET_OBJECTS:${PROJECT_NAME}>
    CACHE INTERNAL "gismo extensions to be included")

# Add Umfpack include directories to G+Smo standard include directories
set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${UMFPACK_INCLUDE_DIR}
  CACHE INTERNAL "gismo include directories")

# Link G+Smo to Umfpack libraries (either dynamically or statically)
set(gismo_LINKER ${gismo_LINKER} ${OpenBLAS_LIBNAME} ${UMFPACK_LIBRARY}
    CACHE INTERNAL "gismo extra linker objects")

# Install gsCoreLibs header files
install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo/gsCoreLibs/
        FILES_MATCHING PATTERN "*.h")

##
## Add dependencies
##
if(BLAS_FOUND)
    # Make gsCoreLibs extension dependent on OpenBLAS library
    add_dependencies(${PROJECT_NAME} openblas)
  endif()
  
if(GISMO_WITH_UMFPACK AND NOT Umfpack_FOUND)
    # Make gsCoreLibs extension dependent on Umfpack library
    add_dependencies(${PROJECT_NAME} umfpack)
endif()
