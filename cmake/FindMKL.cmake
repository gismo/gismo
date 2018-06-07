# - Find Intel MKL
# Find the MKL libraries
#
# Options:
#
#   MKL_STATIC       :   use static linking
#   MKL_MULTI_THREADED:   use multi-threading
#   MKL_SDL           :   Single Dynamic Library interface
#
# This module defines the following variables:
#
#   MKL_FOUND            : True if MKL_INCLUDE_DIR are found
#   MKL_INCLUDE_DIR      : where to find mkl.h, etc.
#   MKL_INCLUDE_DIRS     : set when MKL_INCLUDE_DIR found
#   MKL_LIBRARIES        : the library to link against.
#
# https://software.intel.com/en-us/mkl-linux-developer-guide-selecting-libraries-to-link-with
#
# https://software.intel.com/en-us/mkl-linux-developer-guide-dynamically-selecting-the-interface-and-threading-layer
#

#if (MKL_INCLUDE_DIRS AND MKL_LIBRARIES AND MKL_INTERFACE_LIBRARY AND
#    MKL_SEQUENTIAL_LAYER_LIBRARY AND MKL_CORE_LIBRARY)
#  set (MKL_FIND_QUIETLY TRUE)
#endif()

set(MKL_INTEGER64 ON CACHE BOOL "Enable 64 bit integers on MKL")
set(MKL_MULTI_THREADED ON CACHE BOOL "Enable MKL multithreading")
set(MKL_SDL OFF CACHE BOOL "Enable MKL's  Single Dynamic Library (SDL)")

if(DEFINED INTEL_ROOT AND NOT EXISTS INTEL_ROOT) 
    message(WARNING "The path INTEL_ROOT: ${INTEL_ROOT} does not exist")
endif()

if(NOT DEFINED MKL_ROOT)
  if(DEFINED $ENV{MKL_ROOT})
    set(MKL_ROOT $ENV{MKL_ROOT} CACHE STRING "Folder containing MKL")
  else()
    set(MKL_ROOT ${INTEL_ROOT}/mkl CACHE STRING "Folder containing MKL")
  endif()
  if(NOT EXISTS MKL_ROOT)
    message(WARNING "The path MKL_ROOT: ${MKL_ROOT} does not exist")
  endif()
endif()

include(FindPackageHandleStandardArgs)

# intel64, ia32 # bad for MacOS
if(CMAKE_SIZEOF_VOID_P EQUAL 8) 
    set(MPL_ARCH intel64)
else() 
set(MPL_ARCH intel64)
    set(MPL_ARCH ia32) 
# else: Intel MIC
endif() 

if("x${MPL_ARCH}" STREQUAL "x32")
    SET(ARCH_PREFIX "")
else()
    if(MKL_INTEGER64)
        SET(ARCH_PREFIX "_ilp64")
    else(MKL_INTEGER64)
        SET(ARCH_PREFIX "_lp64")
    endif(MKL_INTEGER64)
endif()

# Find include dir
find_path(MKL_INCLUDE_DIR mkl.h PATHS ${MKL_ROOT}/include)

# Find include directory
#  There is no include folder under linux
if(WIN32)
    find_path(INTEL_INCLUDE_DIR omp.h PATHS ${MKL_ROOT}/include)
    set(MKL_INCLUDE_DIR ${MKL_INCLUDE_DIR} ${INTEL_INCLUDE_DIR})
endif()

# Find libraries

# Handle suffix
if(MKL_STATIC)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_SHARED_LIBRARY_SUFFIX})
endif()

if(MKL_SDL) # Using The Single Dynamic Library (SDL) 
    find_library(MKL_LIBRARY mkl_rt PATHS ${MKL_ROOT}/lib/${MPL_ARCH}/)

else() # MKL is composed by four layers: Interface, Threading, Computational and RTL

    ####################### Computational layer #####################
    find_library(MKL_CORE_LIBRARY mkl_core PATHS ${MKL_ROOT}/lib/${MPL_ARCH}/)
    
    #find_library(MKL_FFT_LIBRARY mkl_cdft_core PATHS ${MKL_ROOT}/lib/${MPL_ARCH}/)
    #find_library(MKL_SCALAPACK_LIBRARY mkl_scalapack${ARCH_PREFIX} PATHS ${MKL_ROOT}/lib/${MPL_ARCH}/)
    #find_library(MKL_BLACS_LIBRARY NAMES mkl_blacs_intelmpi${ARCH_PREFIX} mkl_blacs_openmpi${ARCH_PREFIX} mkl_blacs${ARCH_PREFIX}   PATHS ${MKL_ROOT}/lib/${MPL_ARCH}/)

    #find_library(MPI_LIBRARY mpi_mt PATHS /opt/intel/composer_xe_2013.0.079/mpirt/lib/${MPL_ARCH}/)

    ######################### Interface layer #######################
    find_library(MKL_INTERFACE_LIBRARY mkl_intel${ARCH_PREFIX} mkl_intel_c PATHS ${MKL_ROOT}/lib/${MPL_ARCH}/)

    ######################## Threading layer ########################
    if(MKL_MULTI_THREADED)
      message(STATUS "MKL multi-threading enabled.")
      find_library(MKL_THREADING_LIBRARY mkl_intel_thread PATHS ${MKL_ROOT}/lib/${MPL_ARCH}/)

      #find_library(MKL_THREADING_LIBRARY mkl_pgi_thread
      #find_library(MKL_THREADING_LIBRARY mkl_gnu_thread
      #find_library(MKL_THREADING_LIBRARY mkl_tbb_thread
      #find_library(MKL_THREADING_LIBRARY pthread)
    else()
      find_library(MKL_THREADING_LIBRARY mkl_sequential PATHS ${MKL_ROOT}/lib/${MPL_ARCH}/)
    endif()

    ############################ RTL layer ##########################
    find_library(MKL_RTL_LIBRARY iomp5 iomp5md PATHS ${MKL_ROOT}/lib/${MPL_ARCH} ${INTEL_ROOT}/compiler/lib/${MPL_ARCH} ${MKL_ROOT}/../compiler/lib/${MPL_ARCH})

    message("PATHS ${MKL_ROOT}/lib/${MPL_ARCH}/ ${INTEL_ROOT}/compiler/lib/${MPL_ARCH}/")
    
set(MKL_LIBRARY ${MKL_CORE_LIBRARY} ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_RTL_LIBRARY}
  # ${MKL_FFT_LIBRARY} ${MKL_SCALAPACK_LIBRARY} ${MKL_BLACS_LIBRARY} ${MPI_LIBRARY}
  )
endif()

#message("MKL_CORE_LIBRARY: ${MKL_CORE_LIBRARY}")
#message("MKL_INTERFACE_LIBRARY: ${MKL_INTERFACE_LIBRARY}")
#message("MKL_THREADING_LIBRARY: ${MKL_THREADING_LIBRARY}")
#message("MKL_RTL_LIBRARY: ${MKL_RTL_LIBRARY}")

find_package_handle_standard_args(MKL DEFAULT_MSG MKL_INCLUDE_DIR MKL_LIBRARY)

if(MKL_FOUND)
    message(STATUS "Found Intel MKL libraries.")
    set(MKL_LIBRARIES ${MKL_LIBRARY})
    mark_as_advanced (MKL_LIBRARY) 
  else()
      message(ERROR "THE MKL libraries were not found. Set INTEL_ROOT to the folder containing intel compilers and libraries, or set MKL_ROOT to the folder containing MKL.")
  endif()


#set(MKL_LIBRARY ${MKL_CORE_LIBRARY} ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_RTL_LIBRARY}
#  # ${MKL_FFT_LIBRARY} ${MKL_SCALAPACK_LIBRARY} ${MKL_BLACS_LIBRARY} ${MPI_LIBRARY}

  
