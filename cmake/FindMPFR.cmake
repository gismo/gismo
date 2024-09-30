######################################################################
## FindMPFR.cmake
## This file is part of the G+Smo library. 
##
## Author: Angelos Mantzaflaris 
##
##
## MPFR_FOUND       - system has MPFR lib
## MPFR_INCLUDE_DIR - the MPFR include directory
## MPFR_LIBRARY     - Libraries needed to use MPFR
######################################################################

#unset(MPFR_LIBRARY CACHE)

if (MPFR_INCLUDE_DIR AND MPFR_LIBRARY)
   set(MPFR_FIND_QUIETLY TRUE)
endif ()

set(MPFR_PREFIX "" CACHE PATH "The path to the prefix of an MPFR installation")

find_path(MPFR_INCLUDE_DIR mpfr.h 
	PATHS ${MPFR_PREFIX}/include /usr/include /usr/local/include)

find_library(MPFR_LIBRARY NAMES mpfr libmpfr
	PATHS ${MPFR_PREFIX}/lib /usr/lib /usr/local/lib)

#if(MPFR_INCLUDE_DIR AND MPFR_LIBRARY)
#	get_filename_component(MPFR_LIBRARY_DIR ${MPFR_LIBRARY} PATH)
#   set(MPFR_FOUND TRUE)
#endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPFR DEFAULT_MSG MPFR_INCLUDE_DIR MPFR_LIBRARY)

mark_as_advanced(MPFR_INCLUDE_DIR MPFR_LIBRARY)

if(MPFR_FOUND)
   if(NOT MPFR_FIND_QUIETLY)
      MESSAGE(STATUS "Found MPFR: ${MPFR_LIBRARY}")
   endif()
elseif(MPFR_FOUND)
   if(MPFR_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find MPFR")
   endif()
endif()
