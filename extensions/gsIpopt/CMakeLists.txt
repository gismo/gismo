### CMakeLists.txt ---
##
## Author: Angelos Mantzaflaris
## Copyright (C) 2015 - 2016 RICAM-Linz.
######################################################################

project(gsIpopt)

if(GISMO_WITH_CODIPACK)
  include_directories(${CODIPACK_INCLUDE_DIR})
endif()

#unset(IPOPT_INCLUDE_DIR CACHE)
#unset(IPOPT_LIBRARY     CACHE)

#find_package(BLAS REQUIRED)
find_package(IpOpt QUIET)

set(IPOPT_SOURCES gsOptProblem_.cpp)

if (NOT IPOPT_FOUND)
  set(IPOPT_VER "3.12.6") # 3.12.9 conflicts patching
  message("IpOpt not found (define IpOpt_DIR if it exists on the system).")
  message("Will attempt to download and compile IpOpt ${IPOPT_VER} (do: make IpOpt; cmake .")

  #fetch here!

  include(ExternalProject)
  ExternalProject_Add(IpOpt
    SOURCE_DIR        ${CMAKE_CURRENT_BINARY_DIR}/Ipopt
    CONFIGURE_COMMAND ${CMAKE_CURRENT_BINARY_DIR}/Ipopt/configure --prefix=<INSTALL_DIR> #--libdir=${CMAKE_BINARY_DIR}/lib
    URL http://www.coin-or.org/download/source/Ipopt/Ipopt-${IPOPT_VER}.tgz
    PATCH_COMMAND      patch -t -s -p1 < ${gismo_SOURCE_DIR}/external/ipopt.patch -d ${CMAKE_CURRENT_BINARY_DIR}/Ipopt
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_BINARY_DIR}/Ipopt/ThirdParty/Mumps/mumps.patch ${CMAKE_CURRENT_BINARY_DIR}/Ipopt/mumps.patch
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/Ipopt/ThirdParty/Mumps/get.Mumps
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/Ipopt/ThirdParty/Metis/get.Metis
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/Ipopt/ThirdParty/Blas/get.Blas
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_BINARY_DIR}/lib ${CMAKE_CURRENT_BINARY_DIR}/IpOpt-prefix/lib
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_BINARY_DIR}/lib ${CMAKE_CURRENT_BINARY_DIR}/IpOpt-prefix/lib64
    BUILD_COMMAND ${MAKE}
    UPDATE_COMMAND     ""
    #DOWNLOAD_NO_PROGRESS 0
    )
  set(IPOPT_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/IpOpt-prefix/include/coin CACHE INTERNAL "")
  set(IPOPT_LIBRARY ${CMAKE_BINARY_DIR}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}ipopt${CMAKE_SHARED_LIBRARY_SUFFIX}
    CACHE INTERNAL "")
else()
  #add_custom_target (IpOpt)
endif(NOT IPOPT_FOUND)

include_directories(${IPOPT_INCLUDE_DIR})

add_library(${PROJECT_NAME} OBJECT ${IPOPT_SOURCES})

if (NOT IPOPT_FOUND)
  add_dependencies(${PROJECT_NAME} IpOpt)
endif(NOT IPOPT_FOUND)

set_target_properties(${PROJECT_NAME} PROPERTIES
  COMPILE_DEFINITIONS gismo_EXPORTS
  POSITION_INDEPENDENT_CODE ON
  LINKER_LANGUAGE CXX
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN ON
  FOLDER "G+Smo extensions"  )


# Apply same configuration as G+Smo
include( gsConfig)

set(gismo_EXTENSIONS ${gismo_EXTENSIONS} $<TARGET_OBJECTS:${PROJECT_NAME}>
  CACHE INTERNAL "Gismo extensions to be included")
set(gismo_LINKER ${gismo_LINKER} ${IPOPT_LIBRARY}
  CACHE INTERNAL "${PROJECT_NAME} extra linker objects")

install(DIRECTORY ${PROJECT_SOURCE_DIR}
  DESTINATION include/gismo/gsIpopt/
  FILES_MATCHING PATTERN "*.h")
