######################################################################
## FindIpOpt.cmake
## This file is part of the G+Smo library. 
##
## Author: Angelos Mantzaflaris 
######################################################################

unset(IPOPT_LIBRARY     CACHE)
unset(IPOPT_INCLUDE_DIR CACHE)

find_path(IPOPT_INCLUDE_DIR NAMES IpTNLP.hpp HINTS /usr/include/coin ${CMAKE_CURRENT_BINARY_DIR}/IpOpt-prefix/include/coin ${IpOpt_DIR}/include/coin)

find_library(IPOPT_LIBRARY NAMES ipopt libipopt HINTS ${CMAKE_BINARY_DIR}/lib ${IpOpt_DIR}/lib ${IpOpt_DIR}/lib64 ${IpOpt_DIR}/Ipopt/src/Interfaces/.libs)

if(IPOPT_INCLUDE_DIR AND IPOPT_LIBRARY)
	get_filename_component(IPOPT_LIBRARY_DIR ${IPOPT_LIBRARY} PATH)
   set(IPOPT_FOUND TRUE)
endif()

if(IPOPT_FOUND)
      #message("IPOPT_INCLUDE_DIR: ${IPOPT_INCLUDE_DIR}")
      #message("IPOPT_LIBRARY    : ${IPOPT_LIBRARY}")
   if(NOT IPOPT_FIND_QUIETLY)
      MESSAGE(STATUS "Found IpOpt: ${IPOPT_LIBRARY}")
   endif()
elseif(NOT IPOPT_FOUND)
   set(IpOpt_DIR "IpOpt_DIR-NOTFOUND " CACHE PATH "The path to the IPOPT library")
   if(IPOPT_FIND_REQUIRED)
      message("Could not find IpOpt library.")
      message(FATAL_ERROR "Set the variable IpOpt_DIR and try again.")
   endif()
endif()
