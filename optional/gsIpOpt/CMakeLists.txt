######################################################################
## CMakeLists.txt --- gsIpOpt
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
######################################################################

project(gsIpOpt)

include( gsConfig)

if(GISMO_WITH_CODIPACK)
  include_directories(${CODIPACK_INCLUDE_DIR})
endif()

#unset(IPOPT_INCLUDE_DIR CACHE)
#unset(IPOPT_LIBRARY     CACHE)

#find_package(BLAS REQUIRED)
find_package(IpOpt QUIET)

set(IPOPT_SOURCES gsIpOpt_.cpp)

#if (NOT IPOPT_FOUND)
  set(IPOPT_VER "3.14.10")
  #message("IpOpt not found (define IpOpt_DIR if it exists on the system).")

  #fetch here!
  include(gsFetch)
  gismo_fetch_directory(IpOpt
    URL https://github.com/coin-or/Ipopt/archive/refs/tags/releases/3.14.10.zip
    DESTINATION    external)

  set(IPOPT_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/IpOpt-prefix/include/coin-or CACHE INTERNAL "")
  set(IPOPT_LIBRARY ${CMAKE_CURRENT_BINARY_DIR}/IpOpt-prefix/lib/${CMAKE_SHARED_LIBRARY_PREFIX}ipopt${CMAKE_SHARED_LIBRARY_SUFFIX}
    CACHE INTERNAL "")

  include(ExternalProject)
  ExternalProject_Add(IpOpt
    BINARY_DIR        ${CMAKE_CURRENT_BINARY_DIR}/IpOpt
    SOURCE_DIR        ${gismo_externals}/IpOpt
    CONFIGURE_COMMAND ${gismo_externals}/IpOpt/configure --prefix=<INSTALL_DIR> #ADD_FCFLAGS=-fallow-argument-mismatch
    #URL http://www.coin-or.org/download/source/Ipopt/Ipopt-${IPOPT_VER}.tgz
#    PATCH_COMMAND      patch -t -s -p1 < ${gismo_SOURCE_DIR}/external/ipopt.patch -d ${CMAKE_CURRENT_BINARY_DIR}/Ipopt
#    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_BINARY_DIR}/Ipopt/ThirdParty/Mumps/mumps.patch ${CMAKE_CURRENT_BINARY_DIR}/Ipopt/mumps.patch
#    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/Ipopt/ThirdParty/Mumps/get.Mumps
#    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/Ipopt/ThirdParty/Metis/get.Metis
#    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/Ipopt/ThirdParty/Blas/get.Blas
#    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/Ipopt/ThirdParty/Lapack/get.Lapack
#    BUILD_COMMAND ${MAKE}
    UPDATE_COMMAND     ""
    BUILD_BYPRODUCTS "${IPOPT_LIBRARY}"
    DOWNLOAD_COMMAND  ""
    #DOWNLOAD_NO_PROGRESS 0
    )
#else()
  #add_custom_target (IpOpt)
#endif(NOT IPOPT_FOUND)

## MUMPS (metis?)

set(MUMPS_INST ${CMAKE_CURRENT_BINARY_DIR}/IpOpt-prefix) #install together with IpOpt

  gismo_fetch_directory(ThirdParty-Mumps
    URL https://github.com/coin-or-tools/ThirdParty-Mumps/archive/refs/tags/releases/3.0.7.zip
    DESTINATION    external)

  gismo_fetch_directory(MUMPS
    URL https://files.inria.fr/gismo/gz/MUMPS.tar.gz
    #URL https://github.com/scivision/mumps/archive/refs/heads/main.zip
    DESTINATION    external/ThirdParty-Mumps)

  set(MUMPS_LIBRARY ${MUMPS_INST}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}coinmumps${CMAKE_SHARED_LIBRARY_SUFFIX} )

  # For Mumps: If using GFortran 10, add
  # ADD_FCFLAGS=-fallow-argument-mismatch to compiler flags.
  ExternalProject_Add(ThirdParty-Mumps
      BINARY_DIR        ${CMAKE_CURRENT_BINARY_DIR}/ThirdParty-Mumps
      SOURCE_DIR        ${gismo_externals}/ThirdParty-Mumps
      #COMMAND           ${gismo_externals}/ThirdParty-Mumps/get.Mumps
      CONFIGURE_COMMAND ${gismo_externals}/ThirdParty-Mumps/configure --prefix=${MUMPS_INST} #ADD_FCFLAGS=-fallow-argument-mismatch
      BUILD_BYPRODUCTS "${MUMPS_LIBRARY}"
      DOWNLOAD_COMMAND  ""
      UPDATE_COMMAND    ""
      )

include_directories(${IPOPT_INCLUDE_DIR})
add_library(${PROJECT_NAME} OBJECT ${IPOPT_SOURCES})

#if (NOT IPOPT_FOUND)
  add_dependencies(IpOpt ThirdParty-Mumps)
  add_dependencies(${PROJECT_NAME} IpOpt)
#endif(NOT IPOPT_FOUND)

set_target_properties(${PROJECT_NAME} PROPERTIES
  COMPILE_DEFINITIONS gismo_EXPORTS
  POSITION_INDEPENDENT_CODE ON
  LINKER_LANGUAGE CXX
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN ON
  FOLDER "G+Smo extensions"  )

set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${IPOPT_INCLUDE_DIR}
  CACHE INTERNAL "Gismo include directories")

set(gismo_EXTENSIONS ${gismo_EXTENSIONS} $<TARGET_OBJECTS:${PROJECT_NAME}>
  CACHE INTERNAL "Gismo extensions to be included")
set(gismo_LINKER ${gismo_LINKER} ${IPOPT_LIBRARY} ${MUMPS_LIBRARY}
  CACHE INTERNAL "${PROJECT_NAME} extra linker objects")

#install(DIRECTORY ${PROJECT_SOURCE_DIR}
#  DESTINATION include/gismo/gsIpOpt/
#  FILES_MATCHING PATTERN "*.h")
