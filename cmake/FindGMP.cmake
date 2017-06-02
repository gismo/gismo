######################################################################
## FindGMP.cmake
## This file is part of the G+Smo library. 
##
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - 2015 RICAM-Linz.
##
##
## GMP_FOUND       - system has GMP lib
## GMP_INCLUDE_DIR - the GMP include directory
## GMP_LIBRARY     - Libraries needed to use GMP
## GMPXX_LIBRARY   - Libraries needed to use GMP C++
######################################################################

#unset(GMP_LIBRARY CACHE)
#unset(GMPXX_LIBRARY CACHE)

if (GMP_INCLUDE_DIR AND GMP_LIBRARY AND GMPXX_LIBRARY)
   set(GMP_FIND_QUIETLY TRUE)
endif ()

find_path(GMP_INCLUDE_DIR NAMES gmp.h )
find_library(GMP_LIBRARY   NAMES gmp libgmp) # for static add first libgmp.a
find_library(GMPXX_LIBRARY NAMES gmpxx libgmpxx )
MESSAGE(STATUS "Found GMP: " ${GMP_LIBRARY} " " ${GMPXX_LIBRARY} )

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GMP DEFAULT_MSG GMP_INCLUDE_DIR GMP_LIBRARY)

mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARY)
