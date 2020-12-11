######################################################################
## CMakeLists.txt
## This file is part of the G+Smo library.
## 
## Authors: Menno Deij-van Rijswijk and Angelos Mantzaflaris 
######################################################################

cmake_minimum_required(VERSION 2.8.12)

if(NOT "x${CMAKE_GENERATOR}" STREQUAL "xVisual Studio 14 2015 Win64")
  message(FATAL_ERROR "Rhino plugin: Generator is NOT Visual Studio 14 2015 Win64. Plugin must be disabled.")
endif()

if(NOT ("${CMAKE_BUILD_TYPE}" STREQUAL "Release" OR "${CMAKE_BUILD_TYPE}" STREQUAL "Debug") )
  message(FATAL_ERROR "Rhino plugin: Build type is neither Debug nor Release")
endif()

# Set to false to try out native cmake generators
set(USE_VCXPROJ TRUE)

#message(STATUS "CMAKE_SYSTEM_VERSION ${CMAKE_SYSTEM_VERSION}")
#message(STATUS "CMAKE_VS_WINDOWS_TARGET_PLATFORM_VERSION ${CMAKE_VS_WINDOWS_TARGET_PLATFORM_VERSION}")

project(RhinoTHBSplines)
  
  #set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/plRhino)

  # On first configuration use: -DCMAKE_SYSTEM_VERSION=8.1
  #set(CMAKE_SYSTEM_VERSION 8.1)

  set(CMAKE_MFC_FLAG 2)
  set(RHINO_SDK "C:/Program\ Files/Rhino\ 6.0\ SDK")
  
  #replace /DWIN32 "" in CMAKE_CXX_FLAGS, CMAKE_RC_FLAGS
    
  #set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/plRhino)

if (USE_VCXPROJ)
  
  set(VERSION_MAJOR 1)
  set(VERSION_MINOR 0)
  set(VERSION_BUGFIX 0)

  if (${CMAKE_GENERATOR} MATCHES 2010)
    set (PLATFORM_TOOLSET v100)
  elseif (${CMAKE_GENERATOR} MATCHES 2015)
    set (PLATFORM_TOOLSET v140)
  endif()

  configure_file(Platform.targets.in ${CMAKE_BINARY_DIR}/Platform.targets)
endif()

  include_directories(${RHINO_SDK}/inc)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/GismoSupport)
  
  add_subdirectory(GismoSupport)
  add_subdirectory(GismoExport)
  add_subdirectory(GismoImport)
  add_subdirectory(THBSplines)

