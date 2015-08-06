
set(IpOpt_DIR "IPOPT_DIR-NOTFOUND " CACHE PATH "The path to the IPOPT library")

  find_path(IPOPT_INCLUDE_DIR NAMES IpTNLP.hpp HINTS /usr/include/coin ${CMAKE_CURRENT_BINARY_DIR}/IpOpt-prefix/include/coin ${IPOPT_DIR}/include)

  find_library(IPOPT_LIBRARY NAMES ipopt libipopt HINTS ${CMAKE_BINARY_DIR}/lib ${IPOPT_DIR}/lib)

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
   if(IPOPT_FIND_REQUIRED) ## TO DO, fix
      message("Could not find IpOpt library.")
      message(FATAL_ERROR "Set the variable IpOpt_DIR and try again.")
   endif()
endif()
