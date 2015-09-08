######################################################################
## FindPardiso.cmake ---
## This file is part of the G+Smo library. 
##
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - 2015 RICAM-Linz.
######################################################################


# First try: Pardiso compiled with Intel fortran
if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
   find_library(PARDISO_LIBRARY NAMES pardiso500-INTEL1301-X86-64 pardiso412-INTEL120-X86-64 pardiso411-INTEL101-X86-64 HINTS ${CMAKE_BINARY_DIR}/lib ${Pardiso_DIR} ${Pardiso_DIR}/lib)

if(PARDISO_LIBRARY)

   add_library(Pardiso SHARED IMPORTED)
   set_property(TARGET Pardiso PROPERTY IMPORTED_LOCATION ${PARDISO_LIBRARY})
   #set_property(TARGET Pardiso PROPERTY IMPORTED_IMPLIB   ${PARDISO_LIBRARY.lib}) #ms windows
   set_property(TARGET Pardiso PROPERTY IMPORTED_LINK_INTERFACE_LANGUAGES "Fortran")
   #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lifcore") #-lgfortran 

#   find_package(LAPACK REQUIRED)
#   set_target_properties(Pardiso PROPERTIES IMPORTED_LINK_INTERFACE_LIBRARIES "${LAPACK_LIBRARIES}")

  # Note libmkl_intel.so or libmkl_intel_lp64.so contains "pardiso/pardiso64"

   set(PARDISO_FOUND TRUE)

endif(PARDISO_LIBRARY)

endif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")

# Second try: Pardiso compiled with GNU GCC
if(NOT PARDISO_FOUND)

   find_library(PARDISO_LIBRARY NAMES pardiso500-GNU481-X86-64 pardiso500-GNU472-X86-64 pardiso412-GNU450-X86-64 pardiso412-GNU430-X86-64 pardiso411-GNU443-X86-64 pardiso500-MACOS-X86-64 libpardiso500-WIN-X86-64 libpardiso412-WIN-X86-64 libpardiso412-WIN-X86
             HINTS ${CMAKE_BINARY_DIR}/lib ${Pardiso_DIR} ${Pardiso_DIR}/lib)

  if(PARDISO_LIBRARY)
     add_library(Pardiso SHARED IMPORTED)
     #set_property(TARGET Pardiso PROPERTY IMPORTED_IMPLIB   ${PARDISO_LIBRARY.lib}) #ms windows
     set_property(TARGET Pardiso PROPERTY IMPORTED_LOCATION ${PARDISO_LIBRARY})
     set(PARDISO_FOUND TRUE)

     find_package(LAPACK REQUIRED)
     set_target_properties(Pardiso PROPERTIES IMPORTED_LINK_INTERFACE_LIBRARIES "${LAPACK_LIBRARIES}")

   #find_package(OpenMP)
   #set_target_properties(Pardiso IMPORTED_LINK_INTERFACE_LIBRARIES ${OpenMP_CXX_FLAGS})

  endif(PARDISO_LIBRARY)

endif(NOT PARDISO_FOUND)

# Third try: Pardiso from Intel MKL
if(NOT PARDISO_FOUND)
 
  # see https://software.intel.com/en-us/articles/link-to-intel-mkl-sparse-solvers
  # and https://software.intel.com/en-us/node/470282
  # and https://software.intel.com/en-us/articles/intel-mkl-pardiso
  find_package(MKL QUIET)
  if(MKL_LIBRARIES)
     add_library(Pardiso SHARED IMPORTED)
     set_target_properties(Pardiso PROPERTIES IMPORTED_LINK_INTERFACE_LIBRARIES "${MKL_LIBRARIES}")
     set(PARDISO_FOUND TRUE)
  endif(MKL_LIBRARIES)

endif(NOT PARDISO_FOUND)

if(PARDISO_FOUND)
   get_filename_component(PARDISO_LIBRARY_DIR ${PARDISO_LIBRARY} PATH)
      #message("PARDISO_LIBRARY    : ${PARDISO_LIBRARY}")
   if(NOT PARDISO_FIND_QUIETLY)
      MESSAGE(STATUS "Found Pardiso: ${PARDISO_LIBRARY}")
   endif()
else(PARDISO_FOUND)
   set(Pardiso_DIR "Pardiso_DIR-NOTFOUND " CACHE PATH "The path to the PARDISO library")
   if(Pardiso_FIND_REQUIRED)
      message("Could not find PARDISO library.")
      message(FATAL_ERROR "Set the variable Pardiso_DIR and try again.")
   endif()
endif(PARDISO_FOUND)
