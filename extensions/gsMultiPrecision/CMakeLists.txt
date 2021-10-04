######################################################################
### CMakeLists.txt --- gsMultiPrecision
## This file is part of the G+Smo library.
## 
## Author: Angelos Mantzaflaris, Matthias Moller
######################################################################

# Multi-precision floating-point extension
project(gsMultiPrecisionExtension)

# Use internally built GMP/MPFR libraries
set(GISMO_BUILD_MULTIPRECISION BOOL ON CACHE INTERNAL "")

# Collect file names
aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HEADERS)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_SOURCES)
aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HPPFILES)

# Apply same configuration as G+Smo
include(gsConfig)

##
## GMP library
##
if (GISMO_WITH_MPFR OR GISMO_WITH_GMP)
  # Look for pre-installed GMP library
  find_package(GMP QUIET)
  
  if (NOT GMP_FOUND)
    # Set GMP version
    set(GMP_VER "6.2.1")  
    
    # Download GMP sources at configure time
    include(gsFetch)
    gismo_fetch_directory(gmp
        URL            https://gmplib.org/download/gmp/gmp-${GMP_VER}.tar.bz2
        DESTINATION    external
        )

    # Build GMP library at compile time
    include(ExternalProject)
    ExternalProject_Add(gmp
        BINARY_DIR           ${CMAKE_CURRENT_BINARY_DIR}/gmp
        SOURCE_DIR           ${gismo_externals}/gmp
        CONFIGURE_COMMAND    CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ${gismo_externals}/gmp/configure --enable-cxx --enable-shared=no --enable-static=yes --with-pic --prefix=<INSTALL_DIR>
        DOWNLOAD_COMMAND     ""
        UPDATE_COMMAND       ""
        )
    
    # Set GMP libraries
    set(GMP_LIBRARY ${CMAKE_CURRENT_BINARY_DIR}/gmp-prefix/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gmp${CMAKE_STATIC_LIBRARY_SUFFIX} CACHE INTERNAL "")
    set(GMPXX_LIBRARY ${CMAKE_CURRENT_BINARY_DIR}/gmp-prefix/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gmpxx${CMAKE_STATIC_LIBRARY_SUFFIX} CACHE INTERNAL "")
    
    # Set GMP library and include directories
    ExternalProject_Get_Property(gmp install_dir)
    set(GMP_LIBRARY_DIR ${install_dir}/lib CACHE INTERNAL "")
    set(GMP_INCLUDE_DIR ${install_dir}/include CACHE INTERNAL "")
    include_directories(${GMP_INCLUDE_DIR})
    
    # Install GMP header files
    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gmp-prefix/include
        DESTINATION include/gismo/
        FILES_MATCHING PATTERN "*.h")
  endif (NOT GMP_FOUND)
endif (GISMO_WITH_MPFR OR GISMO_WITH_GMP)

##
## MPFR library
##
if (GISMO_WITH_MPFR)
  # Look for pre-installed MPFR library
  find_package(MPFR QUIET)
  
  if (NOT MPFR_FOUND)  
    # Set MPFR version
    set(MPFR_VER "4.1.0")
    
    # Download MPFR sources at configure time
    include(gsFetch)
    gismo_fetch_directory(mpfr
        URL            https://www.mpfr.org/mpfr-${MPFR_VER}/mpfr-${MPFR_VER}.tar.gz
        DESTINATION    external
        )  
    
    if (NOT GMP_LIBRARY_DIR)
        get_filename_component(GMP_LIBRARY_DIR ${GMP_LIBRARY} PATH)
    endif()
    
    # Build MPFR library at compile time
    include(ExternalProject)
    ExternalProject_Add(mpfr
        BINARY_DIR           ${CMAKE_CURRENT_BINARY_DIR}/mpfr
        SOURCE_DIR           ${gismo_externals}/mpfr
        CONFIGURE_COMMAND    CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ${gismo_externals}/mpfr/configure --with-gmp-include=${GMP_INCLUDE_DIR} --with-gmp-lib=${GMP_LIBRARY_DIR} --enable-shared=no --enable-static=yes  --with-pic --prefix=<INSTALL_DIR>
        DOWNLOAD_COMMAND     ""
        UPDATE_COMMAND       ""
        )
    
    # Set MPFR libraries
    set(MPFR_LIBRARY ${CMAKE_CURRENT_BINARY_DIR}/mpfr-prefix/lib/${CMAKE_STATIC_LIBRARY_PREFIX}mpfr${CMAKE_STATIC_LIBRARY_SUFFIX} CACHE INTERNAL "")
    
    # Set MPFR library and include directories
    ExternalProject_Get_Property(mpfr install_dir)
    set(MPFR_LIBRARY_DIR ${install_dir}/lib CACHE INTERNAL "")
    set(MPFR_INCLUDE_DIR ${install_dir}/include CACHE INTERNAL "")
    include_directories(${MPFR_INCLUDE_DIR})
    
    # Install MPFR header files
    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mpfr-prefix/include
        DESTINATION include/gismo/
        FILES_MATCHING PATTERN "*.h")
  endif(NOT MPFR_FOUND)
endif(GISMO_WITH_MPFR)

##
## gsMultiPrecision extension
##

# Compile gsMultiPrecision extension as part of the G+Smo library
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

# Add gsMultiPrecision extension to the list of G+Smo extensions
set(gismo_EXTENSIONS ${gismo_EXTENSIONS} $<TARGET_OBJECTS:${PROJECT_NAME}>
    CACHE INTERNAL "Gismo extensions to be included")

# Add GMP and MPFR include directories to G+Smo standard include directories
set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${GMP_INCLUDE_DIR} ${MPFR_INCLUDE_DIR}
  CACHE INTERNAL "Gismo include directories")

# Link G+Smo to GMP, GMPXX and MPFR libraries (either dynamically or statically)
set(gismo_LINKER ${gismo_LINKER} ${MPFR_LIBRARY} ${GMPXX_LIBRARY} ${GMP_LIBRARY}
    CACHE INTERNAL "Gismo extra linker objects")

# Install gsMultiPrecision header files
install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo/gsMultiPrecision/
        FILES_MATCHING PATTERN "*.h")

##
## Add dependencies
##
if ((GISMO_WITH_GMP OR GISMO_WITH_MPFR) AND NOT GMP_FOUND)
    # Make gsMultiPrecision extension dependent on GMP library
    add_dependencies(${PROJECT_NAME} gmp)
endif()

if (GISMO_WITH_MPFR AND NOT MPFR_FOUND)
    # Make gsMultiPrecision extension dependent on MPFR library
    add_dependencies(${PROJECT_NAME} mpfr)

    # Make MPFR library dependent on GMP library if it is built from sources
    if (NOT GMP_FOUND)
        add_dependencies(mpfr gmp)
    endif(NOT GMP_FOUND)
endif()
