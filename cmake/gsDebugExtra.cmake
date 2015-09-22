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
    if(GISMO_CLANG_DEBUG)
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
    endif()
  elseif ("x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xGNU")
    # using GCC

    string(REPLACE "-DNDEBUG" "" CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
    string(REPLACE "-DNDEBUG" "" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
    string(REPLACE "-DNDEBUG" "" CMAKE_CXX_FLAGS_RELWITHDEBINFO ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
    string(REPLACE "-DNDEBUG" "" CMAKE_CXX_FLAGS_MINSIZEREL ${CMAKE_CXX_FLAGS_MINSIZEREL})
    string(REPLACE "-DNDEBUG" "" CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG})
    string(REPLACE "-DNDEBUG" "" CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE})
    string(REPLACE "-DNDEBUG" "" CMAKE_C_FLAGS_RELWITHDEBINFO ${CMAKE_C_FLAGS_RELWITHDEBINFO})
    string(REPLACE "-DNDEBUG" "" CMAKE_C_FLAGS_MINSIZEREL ${CMAKE_C_FLAGS_MINSIZEREL})

    #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -rdynamic")
    # Enable checked iterators
    add_definitions(-D_GLIBCXX_DEBUG)
  elseif ("x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xIntel")
    # using Intel C++
  elseif ("x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xMSVC")
    # using Visual Studio C++
    # Enable checked iterators
    string(REPLACE "/D_SECURE_SCL=0" "" CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
    string(REPLACE "/D_SECURE_SCL=0" "" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
    string(REPLACE "/D_SECURE_SCL=0" "" CMAKE_CXX_FLAGS_RELWITHDEBINFO ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
    string(REPLACE "/D_SECURE_SCL=0" "" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_MINSIZEREL})
    FIND_PACKAGE(DbgHelp)

    string(REPLACE "/D NDEBUG" "" CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
    string(REPLACE "/D NDEBUG" "" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
    string(REPLACE "/D NDEBUG" "" CMAKE_CXX_FLAGS_RELWITHDEBINFO ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
    string(REPLACE "/D NDEBUG" "" CMAKE_CXX_FLAGS_MINSIZEREL ${CMAKE_CXX_FLAGS_MINSIZEREL})
    string(REPLACE "/D NDEBUG" "" CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG})
    string(REPLACE "/D NDEBUG" "" CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE})
    string(REPLACE "/D NDEBUG" "" CMAKE_C_FLAGS_RELWITHDEBINFO ${CMAKE_C_FLAGS_RELWITHDEBINFO})
    string(REPLACE "/D NDEBUG" "" CMAKE_C_FLAGS_MINSIZEREL ${CMAKE_C_FLAGS_MINSIZEREL})

  endif()

  #add_definitions(-DEIGEN_INTERNAL_DEBUGGING)

endif(GISMO_EXTRA_DEBUG)
