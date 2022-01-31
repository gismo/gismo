######################################################################
### CMakeLists.txt --- gsMultiPrecision
## This file is part of the G+Smo library.
## 
## Author: Angelos Mantzaflaris, Matthias Moller
######################################################################

# Look for pre-installed GMP library
find_package(GMP) #QUIET

if (NOT GMP_FOUND)
  # Set GMP version
  set(GMP_VER "6.2.1")  
  
  # Download GMP sources at configure time
  include(gsFetch)
  gismo_fetch_directory(gmp
    URL            https://gmplib.org/download/gmp/gmp-${GMP_VER}.tar.bz2
    DESTINATION    external
    )

    # Set GMP libraries
  set(GMP_LIBRARY ${CMAKE_BINARY_DIR}/gmp-prefix/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gmp${CMAKE_STATIC_LIBRARY_SUFFIX} CACHE INTERNAL "")
  set(GMPXX_LIBRARY ${CMAKE_BINARY_DIR}/gmp-prefix/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gmpxx${CMAKE_STATIC_LIBRARY_SUFFIX} CACHE INTERNAL "")

  # Build GMP library at compile time
  include(ExternalProject)
  ExternalProject_Add(gmp
    BINARY_DIR           ${CMAKE_BINARY_DIR}/gmp
    SOURCE_DIR           ${PROJECT_SOURCE_DIR}/external/gmp
    CONFIGURE_COMMAND    CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ${PROJECT_SOURCE_DIR}/external/gmp/configure --enable-cxx --enable-shared=no --enable-static=yes --with-pic --prefix=${CMAKE_BINARY_DIR}/gmp-prefix
    DOWNLOAD_COMMAND     ""
    UPDATE_COMMAND       ""
    BUILD_BYPRODUCTS     "${GMP_LIBRARY};${GMPXX_LIBRARY}"
    )

  # Set GMP library and include directories
  #ExternalProject_Get_Property(gmp install_dir)
  #message("GMP directory: ${install_dir}")
  set(GMP_LIBRARY_DIR "${CMAKE_BINARY_DIR}/gmp-prefix/lib" CACHE INTERNAL "")
  set(GMP_INCLUDE_DIR "/.${CMAKE_BINARY_DIR}/gmp-prefix/include" CACHE INTERNAL "") #note: prefix is to gix a bug with relative paths for ninja
  include_directories(${GMP_INCLUDE_DIR})

  # Install GMP header files
  install(DIRECTORY ${CMAKE_BINARY_DIR}/gmp-prefix/include
    DESTINATION include/gismo/
    FILES_MATCHING PATTERN "*.h")

  #add_dependencies(gismo mpfr)
endif (NOT GMP_FOUND)

# Add GMP and MPFR include directories to G+Smo standard include directories
set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${GMP_INCLUDE_DIR}
  CACHE INTERNAL "gismo include directories" FORCE)

# Link G+Smo to GMP, GMPXX and MPFR libraries (either dynamically or statically)
set(gismo_LINKER ${gismo_LINKER} ${GMPXX_LIBRARY} ${GMP_LIBRARY}
    CACHE INTERNAL "Gismo extra linker objects")

set(GISMO_EXTERNALS ${GISMO_EXTERNALS} "gsGmp"
  CACHE INTERNAL "List of externals" FORCE)

# Install gsGmp header files
#install(DIRECTORY ${PROJECT_SOURCE_DIR}
#        DESTINATION include/gismo/gsGmp/
#        FILES_MATCHING PATTERN "*.h")
