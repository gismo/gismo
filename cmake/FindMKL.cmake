######################################################################
## FindMKL.cmake
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
##
##
## Set INTEL_ROOT to the intel folder in your system
##
## Options:
##
##   MKL_SDL    :  Single Dynamic Library interface
##                 When this option is enabled all the other options have no effect
##   MKL_STATIC :  use static linking
##   MKL_ILP64  :  Use 64bit integers
##
## This module defines the following variables:
##
##   MKL_FOUND            : True if MKL is found
##   MKL_INCLUDE_DIR      : where to find mkl.h
##   MKL_INCLUDE_DIRS     : set when MKL_INCLUDE_DIR found
##   MKL_LIBRARIES        : the library to link against.
##
## Related info:
##
## https://software.intel.com/en-us/mkl-linux-developer-guide-selecting-libraries-to-link-with
## https://software.intel.com/en-us/mkl-linux-developer-guide-dynamically-selecting-the-interface-and-threading-layer
## https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
##
######################################################################

set(MKL_SDL OFF CACHE BOOL "Enable MKL's Single Dynamic Library (SDL)")
set(MKL_STATIC OFF CACHE BOOL "Enable static linkage")
#set(MKL_THREADING CACHE STRING "MKL multithreading backend")
set(MKL_ILP64 OFF CACHE BOOL "Enable 64 bit integers on MKL -- assumes add_definitions(-DMKL_ILP64)") # beware of segfault

unset(MKL_INCLUDE_DIR CACHE)
unset(MKL_LIBRARIES CACHE)
unset(MKL_INCLUDE_DIR CACHE)
unset(MKL_CORE_LIBRARY CACHE)
unset(MKL_INTERFACE_LIBRARY CACHE)
unset(MKL_THREADING_LIBRARY CACHE)
unset(MKL_RTL_LIBRARY CACHE)
unset(MKL_FFT_LIBRARY CACHE)
unset(MKL_SCALAPACK_LIBRARY CACHE)
unset(MKL_BLACS_LIBRARY CACHE)

if(DEFINED INTEL_ROOT AND NOT EXISTS "${INTEL_ROOT}")
    message(WARNING "The path INTEL_ROOT: ${INTEL_ROOT} does not exist")
  else()
    file(GLOB irootL "/*/intel/compilers_and_*" "/*/intel/parallel*/compilers_and_*")
    file(GLOB irootW "C:/Program Files (x86)/IntelSWTools/compilers_and_*" "C:/Program Files (x86)/Intel/Composer*")
    find_path(INTEL_ROOT "mkl" PATHS /opt/intel "${irootL}" "${irootL}/linux" "${irootW}" "${irootW}/windows")
    message(STATUS "Intel root: ${INTEL_ROOT}")
endif()

find_path(MKL_ROOT include/mkl.h PATHS ${INTEL_ROOT}/mkl ${INTEL_ROOT} $ENV{MKLROOT} DOC "Folder contains MKL")

include(FindPackageHandleStandardArgs)

# Find include dir
find_path(MKL_INCLUDE_DIR mkl.h PATHS ${MKL_ROOT} PATH_SUFFIXES include)

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

# intel64, ia32 # bad for macOS
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(MPL_ARCH intel64)
else()
    set(MPL_ARCH ia32)
    # else: Intel MIC
    # set(MPL_ARCH mic)
endif()

if(MKL_SDL) # Using The Single Dynamic Library (SDL)
    find_library(MKL_LIBRARY mkl_rt PATHS ${MKL_ROOT}/lib/${MPL_ARCH})

else() # MKL is composed by four layers: Interface, Threading, Computational and RTL

  if("x${MPL_ARCH}" STREQUAL "x32")
    SET(ARCH_PREFIX "")
  else()
    if(MKL_ILP64)
      SET(ARCH_PREFIX "_ilp64")
    else(MKL_ILP64)
        SET(ARCH_PREFIX "_lp64")
    endif(MKL_ILP64)
  endif()

    ####################### Computational layer #####################
    find_library(MKL_CORE_LIBRARY mkl_core PATHS ${MKL_ROOT}/lib/${MPL_ARCH})
    find_package (Threads) # mkl_core links to pthreads

    ######################### Interface layer #######################
    find_library(MKL_INTERFACE_LIBRARY mkl_intel${ARCH_PREFIX} mkl_intel_c PATHS ${MKL_ROOT}/lib/${MPL_ARCH})

    if(GISMO_WITH_OPENMP)
      ######################## Threading layer ########################
      message(STATUS "MKL multi-threading enabled.")
      find_library(MKL_THREADING_LIBRARY mkl_intel_thread PATHS ${MKL_ROOT}/lib/${MPL_ARCH})
      #find_library(MKL_THREADING_LIBRARY mkl_pgi_thread
      #find_library(MKL_THREADING_LIBRARY mkl_gnu_thread
      #find_library(MKL_THREADING_LIBRARY mkl_tbb_thread

      ############## Runtime-libraries Layer (RTL) ##################
      find_library(MKL_RTL_LIBRARY iomp5 iomp5md PATHS ${INTEL_RTL_ROOT} ${MKL_ROOT}/lib/${MPL_ARCH} ${INTEL_ROOT}/compiler ${INTEL_ROOT}/compiler/lib/${MPL_ARCH} ${MKL_ROOT}/../compiler/lib/${MPL_ARCH} ${INTEL_ROOT}/lib/${MPL_ARCH})

    else(GISMO_WITH_OPENMP)
      message(STATUS "MKL multi-threading is DISABLED (use GISMO_WITH_OPENMP=ON to enable it)")
      find_library(MKL_THREADING_LIBRARY mkl_sequential PATHS ${MKL_ROOT}/lib/${MPL_ARCH})
    endif(GISMO_WITH_OPENMP)

    ####################### Cluster libraries #####################
    if(GISMO_WITH_MPI)
    find_library(MKL_FFT_LIBRARY mkl_cdft_core PATHS ${MKL_ROOT}/lib/${MPL_ARCH})
    find_library(MKL_SCALAPACK_LIBRARY mkl_scalapack${ARCH_PREFIX} PATHS ${MKL_ROOT}/lib/${MPL_ARCH})
    if("x${MKL_BLACS_MPI}" STREQUAL "xINTELMPI")
      find_library(MKL_BLACS_LIBRARY NAMES mkl_blacs_intelmpi${ARCH_PREFIX} PATHS ${MKL_ROOT}/lib/${MPL_ARCH})
    elseif("x${MKL_BLACS_MPI}" STREQUAL "xSGIMPT")
      find_library(MKL_BLACS_LIBRARY NAMES mkl_blacs_sgimpt${ARCH_PREFIX}   PATHS ${MKL_ROOT}/lib/${MPL_ARCH})
    #elseif("x${MKL_BLACS_MPI}" STREQUAL "xMSMPI")
    #elseif("x${MKL_BLACS_MPI}" STREQUAL "xMPICH2")
    else()
      #Note: libmkl_blacs${ARCH_PREFIX} is deprecated!
      find_library(MKL_BLACS_LIBRARY NAMES mkl_blacs_openmpi${ARCH_PREFIX} mkl_blacs_intelmpi${ARCH_PREFIX} PATHS ${MKL_ROOT}/lib/${MPL_ARCH})
    endif()
    else(GISMO_WITH_MPI)
      message(STATUS "MKL cluster parallelization with MPI is DISABLED (use GISMO_WITH_MPI=ON to enable it)")
    endif(GISMO_WITH_MPI)

    # VML: avx512_mic, mc2, mc, acx512, avx, mc3, cmpt, avx2, def
    #find_library(MKL_VML_LIBRARY NAMES mkl_vml_XXX PATHS ${MKL_ROOT}/lib/${MPL_ARCH})
    #find_library(MPI_LIBRARY mpi_mt PATHS /opt/intel/composer_xe_2013.0.079/mpirt/lib/${MPL_ARCH})

set(MKL_LIBRARY ${MKL_CORE_LIBRARY} ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_RTL_LIBRARY}
  ${MKL_FFT_LIBRARY} ${MKL_SCALAPACK_LIBRARY} ${MKL_BLACS_LIBRARY}
  )
endif()

#message("MKL_INCLUDE_DIR ${MKL_INCLUDE_DIR}")
#message("MKL_CORE_LIBRARY: ${MKL_CORE_LIBRARY}")
#message("MKL_INTERFACE_LIBRARY: ${MKL_INTERFACE_LIBRARY}")
#message("MKL_THREADING_LIBRARY: ${MKL_THREADING_LIBRARY}")
#message("MKL_RTL_LIBRARY: ${MKL_RTL_LIBRARY}")
#message("MKL_FFT_LIBRARY: ${MKL_FFT_LIBRARY}")
#message("MKL_SCALAPACK_LIBRARY: ${MKL_SCALAPACK_LIBRARY}")
#message("MKL_BLACS_LIBRARY: ${MKL_BLACS_LIBRARY}")

find_package_handle_standard_args(MKL "The MKL libraries were not found. Set INTEL_ROOT to the folder containing intel compilers and libraries." MKL_INCLUDE_DIR MKL_LIBRARY)

if(MKL_FOUND)
    message(STATUS "Found Intel MKL libraries")
    if(Threads_FOUND)
      set(MKL_LIBRARIES ${MKL_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})
    else()
      set(MKL_LIBRARIES ${MKL_LIBRARY})
    endif()
    mark_as_advanced (MKL_LIBRARY)
endif()
