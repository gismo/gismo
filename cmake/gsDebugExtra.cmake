######################################################################
#
#    Enables extra debugging features.
#
#    This file is part of the G+Smo library.
#
#    This Source Code Form is subject to the terms of the Mozilla Public
#    License, v. 2.0. If a copy of the MPL was not distributed with this
#    file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#    Author(s): Angelos Mantzaflaris
#
######################################################################

if(GISMO_EXTRA_DEBUG)

  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    # using Clang
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # using GCC
    #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -rdynamic")
    # Enable checked iterators
    add_definitions(-D_GLIBCXX_DEBUG)
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    # using Intel C++
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    # using Visual Studio C++
    # Enable checked iterators
    STRING(REPLACE "/D_SECURE_SCL=0" "" CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
    STRING(REPLACE "/D_SECURE_SCL=0" "" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
    FIND_PACKAGE(DbgHelp) 
  endif()
endif(GISMO_EXTRA_DEBUG)
