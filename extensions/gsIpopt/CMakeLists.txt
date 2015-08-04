### CMakeLists.txt ---
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2015 - RICAM-Linz.
######################################################################

## IpOpt extension
# project(gsIpopt) #todo add dep

  find_library(IpOpt_LIB NAMES ipopt libipopt HINTS ${CMAKE_BINARY_DIR}/lib ${IPOPT_DIR}/lib)
  find_path(IpOpt_INCLUDE_DIRS NAMES IpTNLP.hpp HINTS /usr/include/coin ${CMAKE_CURRENT_BINARY_DIR}/IpOpt-prefix/include/coin ${IPOPT_DIR}/include)

  if(IpOpt_LIB)
      #message(WARNING, "Using ipopt system library")
      #find_library(IpOpt_LIB_Mumps NAMES mumps libmumps cmumps libcmumps )
      #find_library(IpOpt_LIB_Metis NAMES metis libmetis )
      #set(IpOpt_INCLUDE_DIRS /usr/include/coin CACHE INTERNAL "IpOpt include path")

  else(IpOpt_LIB)
     #find_library(IpOpt_LIB_Blas NAMES blas libblas cblas atlcblas libatlcblas PATHS /usr/lib64/ NO_DEFAULT_PATH)
     #find_library(IpOpt_LIB_Lapack NAMES lapack liblapack atlclapack libatlclapack )
     #find_package( LAPACK REQUIRED ) # Problematic on ubuntu
     #find_package( BLAS REQUIRED )
     #find_package( DL REQUIRED )
     message(WARNING "Will Download and compile ipopt 3.12.3")
     #message(WARNING, "Found IpOpt_LIB='${IpOpt_LIB}'")
     find_program( patch_cmd "patch" REQUIRED )

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

      set(IpOpt_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/IpOpt-prefix/include/coin
          CACHE INTERNAL "IpOpt include path")

      #add_dependencies(IpOpt mumps metis)
      #add_dependencies(mumps metis)
      #get_property(IpOpt_location TARGET IpOpt PROPERTY LOCATION)
      #add_dependencies(${PROJECT_NAME} IpOpt)
      #set(IpOpt_LIB "${CMAKE_BINARY_DIR}/lib/libipopt.so" CACHE INTERNAL "IpOpt shared library")
      #set(IpOpt_LIB_Mumps "${CMAKE_BINARY_DIR}/lib/libcoinmumps.so" CACHE INTERNAL "IpOpt shared mumps library")
      #set(IpOpt_LIB_Metis "${CMAKE_BINARY_DIR}/lib/libcoinmetis.so" CACHE INTERNAL "IpOpt shared metis library")
    
  endif(IpOpt_LIB)
  #message("IpOpt headers: ${IpOpt_INCLUDE_DIRS}")


# Apply same configuration as G+Smo
#include( gsConfig)

#add_library(gsOpennurbs OBJECT ${OPENNURBS_SOURCES})

#set_target_properties(gsOpennurbs 
#PROPERTIES POSITION_INDEPENDENT_CODE ON
#COMPILE_DEFINITIONS gismo_EXPORTS
#)

#set(gismo_EXTENSIONS ${gismo_EXTENSIONS} $<TARGET_OBJECTS:gsIpopt>
#    CACHE INTERNAL "Gismo extensions to be included")

#install(DIRECTORY ${PROJECT_SOURCE_DIR}
#        DESTINATION include/gismo/gsOpennurbs/
#        FILES_MATCHING PATTERN "*.h")
