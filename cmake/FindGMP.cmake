######################################################################
## FindGMP.cmake
## This file is part of the G+Smo library. 
##
## Author: Angelos Mantzaflaris 
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

set(GMP_PREFIX "" CACHE PATH "The path to the prefix of an GMP installation")

find_path(GMP_INCLUDE_DIR NAMES gmp.h
  PATHS ${GMP_PREFIX}/include /usr/include /usr/local/include)

find_library(GMP_LIBRARY NAMES gmp libgmp
  PATHS ${GMP_PREFIX}/lib /usr/lib /usr/local/lib)

find_library(GMPXX_LIBRARY NAMES gmpxx libgmpxx
  PATHS ${GMP_PREFIX}/lib /usr/lib /usr/local/lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP DEFAULT_MSG GMP_INCLUDE_DIR GMP_LIBRARY GMPXX_LIBRARY)

mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARY GMPXX_LIBRARY)

if(GMP_FOUND)
   if(NOT GMP_FIND_QUIETLY)
      MESSAGE(STATUS "Found GMP: " ${GMP_LIBRARY} " " ${GMPXX_LIBRARY} )
   endif()
elseif(GMP_FOUND)
   if(GMP_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find GMP")
   endif()
endif()
