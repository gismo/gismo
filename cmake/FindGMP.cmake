# Copyright (c) 2008-2010 Kent State University
# Copyright (c) 2011-2012 Texas A&M University
#
# This file is distributed under the MIT License. See the accompanying file
# LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
# and conditions.

# FIXME: How do I find the version of GMP that I want to use?
# What versions are available?

# NOTE: GMP prefix is understood to be the path to the root of the GMP
# installation library.
set(GMP_PREFIX "" CACHE PATH "The path to the prefix of a GMP installation")


find_path(GMP_INCLUDE_DIR gmp.h 
	PATHS ${GMP_PREFIX}/include /usr/include /usr/local/include)

find_library(GMP_LIBRARY NAMES gmp 
	PATHS ${GMP_PREFIX}/lib /usr/lib /usr/local/lib)


if(GMP_INCLUDE_DIR AND GMP_LIBRARY)
	get_filename_component(GMP_LIBRARY_DIR ${GMP_LIBRARY} PATH)
	set(GMP_FOUND TRUE)
endif()

if(GMP_FOUND)
   if(NOT GMP_FIND_QUIETLY)
      MESSAGE(STATUS "Found GMP: ${GMP_LIBRARY}")
   endif()
elseif(GMP_FOUND)
   if(GMP_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find GMP")
   endif()
endif()
