#
# This file is a backport from CMake version 3.7.x. Therefore, 
#

# CMake version 3.1 and below does not set this flag
if (CMAKE_VERSION VERSION_LESS "3.1")
  set(CMAKE_CXX_STANDARD_COMPUTED_DEFAULT ON)
endif()

if ((CMAKE_SYSTEM_NAME STREQUAL "Darwin") AND (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"))

  # Ignore backport if CMake version is more recent
  if (CMAKE_VERSION VERSION_LESS "3.8")
    
    # Apple Clang
    if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.0)
      set(CMAKE_CXX98_STANDARD_COMPILE_OPTION "-std=c++98")
      set(CMAKE_CXX98_EXTENSION_COMPILE_OPTION "-std=gnu++98")
      
      set(CMAKE_CXX11_STANDARD_COMPILE_OPTION "-std=c++11")
      set(CMAKE_CXX11_EXTENSION_COMPILE_OPTION "-std=gnu++11")
    endif()
    
    if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.1)
      set(CMAKE_CXX14_STANDARD_COMPILE_OPTION "-std=c++14")
      set(CMAKE_CXX14_EXTENSION_COMPILE_OPTION "-std=gnu++14")
    elseif(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.1)
      # AppleClang 5.0 knows this flag, but does not set a __cplusplus macro greater than 201103L
      set(CMAKE_CXX14_STANDARD_COMPILE_OPTION "-std=c++1y")
      set(CMAKE_CXX14_EXTENSION_COMPILE_OPTION "-std=gnu++1y")
    endif()
    
    if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.0)
      if (NOT CMAKE_CXX_COMPILER_FORCED)
        if (NOT CMAKE_CXX_STANDARD_COMPUTED_DEFAULT)
          message(FATAL_ERROR "CMAKE_CXX_STANDARD_COMPUTED_DEFAULT should be set for ${CMAKE_CXX_COMPILER_ID} (${CMAKE_CXX_COMPILER}) version ${CMAKE_CXX_COMPILER_VERSION}")
        endif()
        set(CMAKE_CXX_STANDARD_DEFAULT ${CMAKE_CXX_STANDARD_COMPUTED_DEFAULT})
      elseif(NOT DEFINED CMAKE_CXX_STANDARD_DEFAULT)
        # Compiler id was forced so just guess the default standard level.
        set(CMAKE_CXX_STANDARD_DEFAULT 98)
      endif()
    endif()

  endif()
    
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")

  # Ignore backport if CMake version is more recent
  if (CMAKE_VERSION VERSION_LESS "3.8")
    
    # LLVM Clang
    if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 2.1)
      set(CMAKE_CXX98_STANDARD_COMPILE_OPTION "-std=c++98")
      set(CMAKE_CXX98_EXTENSION_COMPILE_OPTION "-std=gnu++98")
    endif()
    
    if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.1)
      set(CMAKE_CXX11_STANDARD_COMPILE_OPTION "-std=c++11")
      set(CMAKE_CXX11_EXTENSION_COMPILE_OPTION "-std=gnu++11")
    elseif(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 2.1)
      set(CMAKE_CXX11_STANDARD_COMPILE_OPTION "-std=c++0x")
      set(CMAKE_CXX11_EXTENSION_COMPILE_OPTION "-std=gnu++0x")
    endif()
    
    if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.5)
      set(CMAKE_CXX14_STANDARD_COMPILE_OPTION "-std=c++14")
      set(CMAKE_CXX14_EXTENSION_COMPILE_OPTION "-std=gnu++14")
    elseif(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.4)
      set(CMAKE_CXX14_STANDARD_COMPILE_OPTION "-std=c++1y")
      set(CMAKE_CXX14_EXTENSION_COMPILE_OPTION "-std=gnu++1y")
    endif()
    
    if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.4)
      if (NOT CMAKE_CXX_COMPILER_FORCED)
        if (NOT CMAKE_CXX_STANDARD_COMPUTED_DEFAULT)
          message(FATAL_ERROR "CMAKE_CXX_STANDARD_COMPUTED_DEFAULT should be set for ${CMAKE_CXX_COMPILER_ID} (${CMAKE_CXX_COMPILER}) version ${CMAKE_CXX_COMPILER_VERSION}")
        endif()
        set(CMAKE_CXX_STANDARD_DEFAULT ${CMAKE_CXX_STANDARD_COMPUTED_DEFAULT})
      elseif(NOT DEFINED CMAKE_CXX_STANDARD_DEFAULT)
        # Compiler id was forced so just guess the default standard level.
        set(CMAKE_CXX_STANDARD_DEFAULT 98)
      endif()
    endif()

  endif()
  
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

  # Ignore backport if CMake version is more recent
  if (CMAKE_VERSION VERSION_LESS "3.8")
    
    # GNU Compiler Collection
    if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.4)
      # Supported since 4.3
      set(CMAKE_CXX98_STANDARD_COMPILE_OPTION "-std=c++98")
      set(CMAKE_CXX98_EXTENSION_COMPILE_OPTION "-std=gnu++98")
    endif()
    
    if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.7)
      set(CMAKE_CXX11_STANDARD_COMPILE_OPTION "-std=c++11")
      set(CMAKE_CXX11_EXTENSION_COMPILE_OPTION "-std=gnu++11")
    elseif (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.4)
      # 4.3 supports 0x variants
      set(CMAKE_CXX11_STANDARD_COMPILE_OPTION "-std=c++0x")
      set(CMAKE_CXX11_EXTENSION_COMPILE_OPTION "-std=gnu++0x")
    endif()
    
    if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9)
      set(CMAKE_CXX14_STANDARD_COMPILE_OPTION "-std=c++14")
      set(CMAKE_CXX14_EXTENSION_COMPILE_OPTION "-std=gnu++14")
    elseif (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8)
      set(CMAKE_CXX14_STANDARD_COMPILE_OPTION "-std=c++1y")
      set(CMAKE_CXX14_EXTENSION_COMPILE_OPTION "-std=gnu++1y")
    endif()
    
    if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.4)
      if (NOT CMAKE_CXX_COMPILER_FORCED)
        if (NOT CMAKE_CXX_STANDARD_COMPUTED_DEFAULT)
          message(FATAL_ERROR "CMAKE_CXX_STANDARD_COMPUTED_DEFAULT should be set for ${CMAKE_CXX_COMPILER_ID} (${CMAKE_CXX_COMPILER}) version ${CMAKE_CXX_COMPILER_VERSION}")
        endif()
        set(CMAKE_CXX_STANDARD_DEFAULT ${CMAKE_CXX_STANDARD_COMPUTED_DEFAULT})
      elseif(NOT DEFINED CMAKE_CXX_STANDARD_DEFAULT)
        # Compiler id was forced so just guess the default standard level.
        set(CMAKE_CXX_STANDARD_DEFAULT 98)
      endif()
    endif()

  endif()
    
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")

  # Ignore backport if CMake version is more recent
  if (CMAKE_VERSION VERSION_LESS "3.8")
    
    # Intel compiler 
    if("x${CMAKE_CXX_SIMULATE_ID}" STREQUAL "xMSVC")
      set(_std -Qstd)
      set(_ext c++)
    else()
      set(_std -std)
      set(_ext gnu++)
    endif()
    
    if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15.0.2)
      set(CMAKE_CXX14_STANDARD_COMPILE_OPTION "${_std}=c++14")
      # todo: there is no gnu++14 value supported; figure out what to do
      set(CMAKE_CXX14_EXTENSION_COMPILE_OPTION "${_std}=c++14")
    elseif (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15.0.0)
      set(CMAKE_CXX14_STANDARD_COMPILE_OPTION "${_std}=c++1y")
      # todo: there is no gnu++14 value supported; figure out what to do
      set(CMAKE_CXX14_EXTENSION_COMPILE_OPTION "${_std}=c++1y")
    endif()
    
    if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 13.0)
      set(CMAKE_CXX11_STANDARD_COMPILE_OPTION "${_std}=c++11")
      set(CMAKE_CXX11_EXTENSION_COMPILE_OPTION "${_std}=${_ext}11")
    elseif (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 12.1)
      set(CMAKE_CXX11_STANDARD_COMPILE_OPTION "${_std}=c++0x")
      set(CMAKE_CXX11_EXTENSION_COMPILE_OPTION "${_std}=${_ext}0x")
    endif()
    
    if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 12.1)
      set(CMAKE_CXX98_STANDARD_COMPILE_OPTION "${_std}=c++98")
      set(CMAKE_CXX98_EXTENSION_COMPILE_OPTION "${_std}=${_ext}98")
    endif()
    
    if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 12.1)
      if (NOT CMAKE_CXX_COMPILER_FORCED)
        if (NOT CMAKE_CXX_STANDARD_COMPUTED_DEFAULT)
          message(FATAL_ERROR "CMAKE_CXX_STANDARD_COMPUTED_DEFAULT should be set for ${CMAKE_CXX_COMPILER_ID} (${CMAKE_CXX_COMPILER}) version ${CMAKE_CXX_COMPILER_VERSION}")
        else()
          set(CMAKE_CXX_STANDARD_DEFAULT ${CMAKE_CXX_STANDARD_COMPUTED_DEFAULT})
        endif()
      elseif (NOT CMAKE_CXX_STANDARD_COMPUTED_DEFAULT)
        # Compiler id was forced so just guess the default standard level.
        set(CMAKE_CXX_STANDARD_DEFAULT 98)
      endif()
    endif()
    
    unset(_std)
    unset(_ext)

  endif()
  
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "PGI")

  # CMake does not yet provide flags for the Portland Group compiler
  
  # The Portland Group
  if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15.0.0)
    set(CMAKE_CXX98_STANDARD_COMPILE_OPTION "$std=c++98")
    set(CMAKE_CXX98_EXTENSION_COMPILE_OPTION "$std=c++98")
  endif()

  if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15.0.0)
    set(CMAKE_CXX11_STANDARD_COMPILE_OPTION "-std=c++11")
    set(CMAKE_CXX11_EXTENSION_COMPILE_OPTION "-std=c++11")
  else()
    set(CMAKE_CXX11_STANDARD_COMPILE_OPTION "-std=c++0x")
    set(CMAKE_CXX11_EXTENSION_COMPILE_OPTION "-std=c++0x")
  endif()
  
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "SunPro")

  # Ignore backport if CMake version is more recent
  if (CMAKE_VERSION VERSION_LESS "3.8")
    
    # Oracle Solaris Studio
    if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.13)
      set(CMAKE_CXX11_STANDARD_COMPILE_OPTION "-std=c++11")
      set(CMAKE_CXX11_EXTENSION_COMPILE_OPTION "-std=c++11")
    endif()
    
    if(NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.13)
      if (NOT CMAKE_CXX_COMPILER_FORCED)
        if (NOT CMAKE_CXX_STANDARD_COMPUTED_DEFAULT)
          message(FATAL_ERROR "CMAKE_CXX_STANDARD_COMPUTED_DEFAULT should be set for ${CMAKE_CXX_COMPILER_ID} (${CMAKE_CXX_COMPILER}) version ${CMAKE_CXX_COMPILER_VERSION}")
        endif()
        set(CMAKE_CXX_STANDARD_DEFAULT ${CMAKE_CXX_STANDARD_COMPUTED_DEFAULT})
      elseif(NOT DEFINED CMAKE_CXX_STANDARD_DEFAULT)
        # Compiler id was forced so just guess the default standard level.
        set(CMAKE_CXX_STANDARD_DEFAULT 98)
      endif()
    endif()
    
  endif()

else() # MSVC or others
     set(CMAKE_CXX${CMAKE_CXX_STANDARD}_STANDARD_COMPILE_OPTION "")
     set(CMAKE_CXX${CMAKE_CXX_STANDARD}_EXTENSION_COMPILE_OPTION "")     
endif()
