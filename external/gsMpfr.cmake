######################################################################
### CMakeLists.txt --- gsMultiPrecision
## This file is part of the G+Smo library.
## 
## Author: Angelos Mantzaflaris, Matthias Moller
######################################################################

# Look for pre-installed MPFR library
find_package(MPFR) #QUIET

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

  # Set MPFR libraries
  set(MPFR_LIBRARY ${CMAKE_BINARY_DIR}/mpfr-prefix/lib/${CMAKE_STATIC_LIBRARY_PREFIX}mpfr${CMAKE_STATIC_LIBRARY_SUFFIX} CACHE INTERNAL "")

  # Build MPFR library at compile time
  include(ExternalProject)
  ExternalProject_Add(mpfr
    BINARY_DIR           ${CMAKE_BINARY_DIR}/mpfr
    SOURCE_DIR           ${PROJECT_SOURCE_DIR}/external/mpfr
    CONFIGURE_COMMAND    CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ${PROJECT_SOURCE_DIR}/external/mpfr/configure --with-gmp-include=${GMP_INCLUDE_DIR} --with-gmp-lib=${GMP_LIBRARY_DIR} --enable-shared=no --enable-static=yes  --with-pic --prefix=${CMAKE_BINARY_DIR}/mpfr-prefix
    DOWNLOAD_COMMAND     ""
    UPDATE_COMMAND       ""
    BUILD_BYPRODUCTS     "${MPFR_LIBRARY}"
    )

  # Set MPFR library and include directories
  #ExternalProject_Get_Property(mpfr install_dir)
  #message("MPFR directory: ${install_dir}")
  set(MPFR_LIBRARY_DIR ${CMAKE_BINARY_DIR}/mpfr-prefix/lib CACHE INTERNAL "")
  set(MPFR_INCLUDE_DIR "/.${CMAKE_BINARY_DIR}/mpfr-prefix/include" CACHE INTERNAL "")
  include_directories(${MPFR_INCLUDE_DIR})
  
  # Install MPFR header files
  install(DIRECTORY ${CMAKE_BINARY_DIR}/mpfr-prefix/include
    DESTINATION include/gismo/
    FILES_MATCHING PATTERN "*.h")

#  add_dependencies(mpfr gmp)
endif(NOT MPFR_FOUND)

# Add GMP and MPFR include directories to G+Smo standard include directories
set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${MPFR_INCLUDE_DIR}
  CACHE INTERNAL "gismo include directories" FORCE)

# Link G+Smo to GMP, GMPXX and MPFR libraries (either dynamically or statically)
set(gismo_LINKER ${gismo_LINKER} ${MPFR_LIBRARY}
    CACHE INTERNAL "Gismo extra linker objects")

set(GISMO_EXTERNALS ${GISMO_EXTERNALS} "gsMpfr"
  CACHE INTERNAL "List of externals" FORCE)

# Install gsMpfr header files
#install(DIRECTORY ${PROJECT_SOURCE_DIR}
#        DESTINATION include/gismo/gsMpfr/
#        FILES_MATCHING PATTERN "*.h")
