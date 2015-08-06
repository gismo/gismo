### CMakeLists.txt ---
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2015 - RICAM-Linz.
######################################################################

## IpOpt extension
# project(gsIpopt) #todo add dep

find_package(IpOpt QUIET)

if (NOT IPOPT_FOUND)
message("IpOpt not found (define IpOpt_DIR if it exists on the system). Will attempt to download and compile ipopt 3.12.3")
include(ExternalProject)
ExternalProject_Add(IpOpt
          SOURCE_DIR        ${CMAKE_CURRENT_BINARY_DIR}/Ipopt
          CONFIGURE_COMMAND ${CMAKE_CURRENT_BINARY_DIR}/Ipopt/configure --prefix=<INSTALL_DIR> #--libdir=${CMAKE_BINARY_DIR}/lib
          URL http://www.coin-or.org/download/source/Ipopt/Ipopt-3.12.3.tgz
          PATCH_COMMAND      patch -t -s -p1 < ${gismo_SOURCE_DIR}/external/ipopt.patch -d ${CMAKE_CURRENT_BINARY_DIR}/Ipopt
          COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_BINARY_DIR}/Ipopt/ThirdParty/Mumps/mumps.patch ${CMAKE_CURRENT_BINARY_DIR}/Ipopt/mumps.patch
          COMMAND ${CMAKE_CURRENT_BINARY_DIR}/Ipopt/ThirdParty/Mumps/get.Mumps
          COMMAND ${CMAKE_CURRENT_BINARY_DIR}/Ipopt/ThirdParty/Metis/get.Metis
          COMMAND ${CMAKE_CURRENT_BINARY_DIR}/Ipopt/ThirdParty/Blas/get.Blas
          COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_BINARY_DIR}/lib ${CMAKE_CURRENT_BINARY_DIR}/IpOpt-prefix/lib
          COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_BINARY_DIR}/lib ${CMAKE_CURRENT_BINARY_DIR}/IpOpt-prefix/lib64
          BUILD_COMMAND ${MAKE} 
          UPDATE_COMMAND     ""
)

      unset(IPOPT_INCLUDE_DIR CACHE)
      unset(IPOPT_LIBRARY     CACHE)
      set(IPOPT_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/IpOpt-prefix/include/coin)
      set(IPOPT_LIBRARY libipopt)

else()
    add_custom_target (IpOpt)
endif(NOT IPOPT_FOUND)


#message("IPOPT_LIBRARY ${IPOPT_LIBRARY}")
#message("IPOPT_INCLUDE_DIR ${IPOPT_INCLUDE_DIR}")
