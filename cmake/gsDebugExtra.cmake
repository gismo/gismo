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

  if ("x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xClang")
    # using Clang
  elseif ("x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xGNU")
    # using GCC
    #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -rdynamic")
    # Enable checked iterators
    add_definitions(-D_GLIBCXX_DEBUG)
  elseif ("x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xIntel")
    # using Intel C++
  elseif ("x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xMSVC")
    # using Visual Studio C++
    # Enable checked iterators
    STRING(REPLACE "/D_SECURE_SCL=0" "" CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
    STRING(REPLACE "/D_SECURE_SCL=0" "" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
    FIND_PACKAGE(DbgHelp) 
  endif()
endif(GISMO_EXTRA_DEBUG)
