# Copyright (c) 2008-2010 Kent State University
# Copyright (c) 2011-2012 Texas A&M University
#
# This file is distributed under the MIT License. See the accompanying file
# LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
# and conditions.

# FIXME: How do I find the version of MPFR that I want to use?
# What versions are available?

# NOTE: MPFR prefix is understood to be the path to the root of the MPFR
# installation library.
set(MPFR_PREFIX "" CACHE PATH "The path to the previx of an MPFR installation")

find_path(MPFR_INCLUDE_DIR mpfr.h 
	PATHS ${MPFR_PREFIX}/include /usr/include /usr/local/include)

find_library(MPFR_LIBRARY NAMES mpfr 
	PATHS ${MPFR_PREFIX}/lib /usr/lib /usr/local/lib)

if(MPFR_INCLUDE_DIR AND MPFR_LIBRARY)
	get_filename_component(MPFR_LIBRARY_DIR ${MPFR_LIBRARY} PATH)
   set(MPFR_FOUND TRUE)
endif()


if(MPFR_FOUND)
   if(NOT MPFR_FIND_QUIETLY)
      MESSAGE(STATUS "Found MPFR: ${MPFR_LIBRARY}")
   endif()
elseif(MPFR_FOUND)
   if(MPFR_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find MPFR")
   endif()
endif()
