######################################################################
## AddCXXConpileOptions.cmake
## This file is part of the G+Smo library. 
##
## Authors: M. Moeller and A. Mantzaflaris 
######################################################################

set(CMAKE_CXX_STANDARD_DEFAULT 11)

if (CMAKE_CXX_COMPILER_ID STREQUAL "PGI")

  # CMake does not yet provide flags for the Portland Group compiler
  
  # The Portland Group
  if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15.0.0)
    set(CMAKE_CXX98_STANDARD_COMPILE_OPTION "$std=c++98")
    set(CMAKE_CXX98_EXTENSION_COMPILE_OPTION "$std=c++98")
    set(CMAKE_CXX11_STANDARD_COMPILE_OPTION "-std=c++11")
    set(CMAKE_CXX11_EXTENSION_COMPILE_OPTION "-std=c++11")
    set(CMAKE_CXX_STANDARD_DEFAULT 11)
  else()
    set(CMAKE_CXX11_STANDARD_COMPILE_OPTION "-std=c++0x")
    set(CMAKE_CXX11_EXTENSION_COMPILE_OPTION "-std=c++0x")
    set(CMAKE_CXX_STANDARD_DEFAULT 98)
  endif()
  
endif()

if (CMAKE_VERSION VERSION_LESS "3.1")

if ((CMAKE_SYSTEM_NAME STREQUAL "Darwin") AND (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"))

    #also: -stdlib=libc++ 
    
    # Apple Clang
    if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.0)
      set(CMAKE_CXX98_STANDARD_COMPILE_OPTION "-std=c++98")
      set(CMAKE_CXX98_EXTENSION_COMPILE_OPTION "-std=gnu++98")      
      set(CMAKE_CXX11_STANDARD_COMPILE_OPTION "-std=c++11")
      set(CMAKE_CXX11_EXTENSION_COMPILE_OPTION "-std=gnu++11")
      set(CMAKE_CXX_STANDARD_DEFAULT 11)
    endif()
    
    if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.1)
      set(CMAKE_CXX14_STANDARD_COMPILE_OPTION "-std=c++14")
      set(CMAKE_CXX14_EXTENSION_COMPILE_OPTION "-std=gnu++14")
      set(CMAKE_CXX_STANDARD_DEFAULT 14)
    elseif(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.1)
      # AppleClang 5.0 knows this flag, but does not set a __cplusplus macro greater than 201103L
      set(CMAKE_CXX14_STANDARD_COMPILE_OPTION "-std=c++1y")
      set(CMAKE_CXX14_EXTENSION_COMPILE_OPTION "-std=gnu++1y")
      set(CMAKE_CXX_STANDARD_DEFAULT 14)
    endif()
        
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")

    # LLVM Clang
    if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 2.1)
      set(CMAKE_CXX98_STANDARD_COMPILE_OPTION "")
      set(CMAKE_CXX98_EXTENSION_COMPILE_OPTION "-std=gnu++98")
      set(CMAKE_CXX_STANDARD_DEFAULT 98)
    endif()
    
    if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.1)
      set(CMAKE_CXX11_STANDARD_COMPILE_OPTION "-std=c++11")
      set(CMAKE_CXX11_EXTENSION_COMPILE_OPTION "-std=gnu++11")
      #set(CMAKE_CXX_STANDARD_DEFAULT 11) # travis/gcc4.6
    elseif(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 2.1)
      set(CMAKE_CXX11_STANDARD_COMPILE_OPTION "-std=c++0x")
      set(CMAKE_CXX11_EXTENSION_COMPILE_OPTION "-std=gnu++0x")
      #set(CMAKE_CXX_STANDARD_DEFAULT 11)
    endif()
    
    if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.5)
      set(CMAKE_CXX14_STANDARD_COMPILE_OPTION "-std=c++14")
      set(CMAKE_CXX14_EXTENSION_COMPILE_OPTION "-std=gnu++14")
      set(CMAKE_CXX_STANDARD_DEFAULT 14)
    elseif(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.4)
      set(CMAKE_CXX14_STANDARD_COMPILE_OPTION "-std=c++1y")
      set(CMAKE_CXX14_EXTENSION_COMPILE_OPTION "-std=gnu++1y")
      # .. additionally requires gnu libstdc++  greater than 4.6
      # set(CMAKE_CXX_STANDARD_DEFAULT 14)
    endif()

elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

    # GNU Compiler Collection
    if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.4)
      # Supported since 4.3
      set(CMAKE_CXX98_STANDARD_COMPILE_OPTION "-std=c++98")
      set(CMAKE_CXX98_EXTENSION_COMPILE_OPTION "-std=gnu++98")
      set(CMAKE_CXX_STANDARD_DEFAULT 98)
    endif()
    
    if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.7)
      set(CMAKE_CXX11_STANDARD_COMPILE_OPTION "-std=c++11")
      set(CMAKE_CXX11_EXTENSION_COMPILE_OPTION "-std=gnu++11")
      set(CMAKE_CXX_STANDARD_DEFAULT 11)
    elseif (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.4)
      # 4.3 supports 0x variants
      set(CMAKE_CXX11_STANDARD_COMPILE_OPTION "-std=c++0x")
      set(CMAKE_CXX11_EXTENSION_COMPILE_OPTION "-std=gnu++0x")
      set(CMAKE_CXX_STANDARD_DEFAULT 14)
    endif()
    
    if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9)
      set(CMAKE_CXX14_STANDARD_COMPILE_OPTION "-std=c++14")
      set(CMAKE_CXX14_EXTENSION_COMPILE_OPTION "-std=gnu++14")
      set(CMAKE_CXX_STANDARD_DEFAULT 14)
    elseif (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8)
      set(CMAKE_CXX14_STANDARD_COMPILE_OPTION "-std=c++1y")
      set(CMAKE_CXX14_EXTENSION_COMPILE_OPTION "-std=gnu++1y")
      set(CMAKE_CXX_STANDARD_DEFAULT 14)
    endif()
   
elseif ( "x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xIntel")

    # Intel compiler 
    if("x${CMAKE_CXX_SIMULATE_ID}" STREQUAL "xMSVC")
      set(_std -Qstd)
      set(_ext c++)
    else()
      set(_std -std)
      set(_ext gnu++)
    endif()

    if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 13.1)
      set(CMAKE_CXX98_STANDARD_COMPILE_OPTION "${_std}=c++98")
      set(CMAKE_CXX98_EXTENSION_COMPILE_OPTION "${_std}=${_ext}98")
      set(CMAKE_CXX_STANDARD_DEFAULT 98)
    endif()

    if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15.0.2)
      set(CMAKE_CXX11_STANDARD_COMPILE_OPTION "${_std}=c++11")
      set(CMAKE_CXX11_EXTENSION_COMPILE_OPTION "${_std}=${_ext}11")
      set(CMAKE_CXX14_STANDARD_COMPILE_OPTION "${_std}=c++14")
      # todo: there is no gnu++14 value supported; figure out what to do
      set(CMAKE_CXX14_EXTENSION_COMPILE_OPTION "${_std}=c++14")
      set(CMAKE_CXX_STANDARD_DEFAULT 14)
    elseif (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15.0.0)
      set(CMAKE_CXX11_STANDARD_COMPILE_OPTION "${_std}=c++0x")
      set(CMAKE_CXX11_EXTENSION_COMPILE_OPTION "${_std}=${_ext}0x")
      set(CMAKE_CXX14_STANDARD_COMPILE_OPTION "${_std}=c++1y")
      # todo: there is no gnu++14 value supported; figure out what to do
      set(CMAKE_CXX14_EXTENSION_COMPILE_OPTION "${_std}=c++1y")
      set(CMAKE_CXX_STANDARD_DEFAULT 14)
    endif()
           
    unset(_std)
    unset(_ext)
    
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "SunPro")
    
    # Oracle Solaris Studio
    if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.13)
      set(CMAKE_CXX11_STANDARD_COMPILE_OPTION "-std=c++11")
      set(CMAKE_CXX11_EXTENSION_COMPILE_OPTION "-std=c++11")
      set(CMAKE_CXX_STANDARD_DEFAULT 11)
    endif()

endif()

endif() # cmake 3.1

# Set C++ standard
if(DEFINED GISMO_BUILD_CPP11) # B.C.
  if(${GISMO_BUILD_CPP11}) # B.C.
    set(CMAKE_CXX_STANDARD 11 CACHE INTERNAL "")
  else()
    set(CMAKE_CXX_STANDARD 98 CACHE INTERNAL "")
  endif()
  unset(GISMO_BUILD_CPP11 CACHE)
endif()

if (NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD ${CMAKE_CXX_STANDARD_DEFAULT} CACHE INTERNAL "")
endif()

# Apply for Cmake less than 3.1
if (CMAKE_VERSION VERSION_LESS "3.1")

  if ( NOT "x${CMAKE_CXX_STANDARD}" STREQUAL "x98" AND
       ${CMAKE_CXX_STANDARD_DEFAULT} LESS ${CMAKE_CXX_STANDARD})
      message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} supports at most C++${CMAKE_CXX_STANDARD_DEFAULT} (requested: ${CMAKE_CXX_STANDARD}).")
  endif()

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX${CMAKE_CXX_STANDARD}_STANDARD_COMPILE_OPTION}")
endif()#cmake<3.1


# Bugfix for windows/msvc systems
if(NOT DEFINED CMAKE_CXX${CMAKE_CXX_STANDARD}_STANDARD_COMPILE_OPTION)
      set(CMAKE_CXX${CMAKE_CXX_STANDARD}_STANDARD_COMPILE_OPTION  "")
      set(CMAKE_CXX${CMAKE_CXX_STANDARD}_EXTENSION_COMPILE_OPTION "")
endif()


# if(NOT DEFINED CMAKE_CXX${CMAKE_CXX_STANDARD}_STANDARD_COMPILE_OPTION) 
#   if((CMAKE_SYSTEM_NAME STREQUAL "Darwin") AND (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"))
#     # set(CMAKE_CXX_FLAGS "-stdlib=libc++")
#     set(CMAKE_CXX${CMAKE_CXX_STANDARD}_STANDARD_COMPILE_OPTION "-std=c++${CMAKE_CXX_STANDARD}")
#     set(CMAKE_CXX${CMAKE_CXX_STANDARD}_EXTENSION_COMPILE_OPTION "-std=gnu++${CMAKE_CXX_STANDARD}")
#   elseif((CMAKE_SYSTEM_NAME STREQUAL "Windows") AND MSVC)
#     set(CMAKE_CXX${CMAKE_CXX_STANDARD}_STANDARD_COMPILE_OPTION "")
#     set(CMAKE_CXX${CMAKE_CXX_STANDARD}_EXTENSION_COMPILE_OPTION "")
#   endif()
# endif()

# if (CMAKE_VERSION VERSION_LESS "3.1" AND CMAKE_COMPILER_IS_GNUCC)
#   if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.1)
#     CHECK_CXX_COMPILER_FLAG("-std=c++${CMAKE_CXX_STANDARD}"
#       COMPILER_SUPPORTS_CXX${CMAKE_CXX_STANDARD})
#     if(COMPILER_SUPPORTS_CXX${CMAKE_CXX_STANDARD})
#       if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++${CMAKE_CXX_STANDARD} -stdlib=libc++")
#       else()
#         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++${CMAKE_CXX_STANDARD}")
#       endif()
#     else()
#       CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
#       if(COMPILER_SUPPORTS_CXX0X)
#         if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#           set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -stdlib=libc++")
#         else()
#           set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
#         endif()
#       else()
#         message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++${CMAKE_CXX_STANDARD}.")
#       endif()
#     endif()
#   else() #gcc 6.1
#     if(CMAKE_CXX_STANDARD EQUAL "98")
#       set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++98")
#     endif() 
#   endif() 
# endif()#cmake<3.1

