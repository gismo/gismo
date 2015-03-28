# - Find Intel MKL
# Find the MKL libraries
#
# Options:
#
#   MKL_STATAIC       :   use static linking
#   MKL_MULTI_THREADED:   use multi-threading
#   MKL_SDL           :   Single Dynamic Library interface
#
# This module defines the following variables:
#
#   MKL_FOUND            : True if MKL_INCLUDE_DIR are found
#   MKL_INCLUDE_DIR      : where to find mkl.h, etc.
#   MKL_INCLUDE_DIRS     : set when MKL_INCLUDE_DIR found
#   MKL_LIBRARIES        : the library to link against.


#message("CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES}")


include(FindPackageHandleStandardArgs)

set(INTEL_ROOT "/opt/intel" CACHE PATH "Folder contains intel libs")
set(MKL_ROOT ${INTEL_ROOT}/mkl CACHE PATH "Folder contains MKL")

# intel64, ia32 # bad for MacOS
if(CMAKE_SIZEOF_VOID_P EQUAL 8) 
    set(MPL_ARCH intel64)
else() 
set(MPL_ARCH intel64)
    set(MPL_ARCH ia32) 
endif() 

#set(MKL_INTEGER64 OFF CACHE BOOL "Enable 64 bit integers on MKL")
set(MKL_INTEGER64 OFF)

if(variable STREQUAL "32")
    SET(ARCH_PREFIX "")
else(variable STREQUAL "32")
    if(MKL_INTEGER64)
        SET(ARCH_PREFIX "_ilp64")
    else(MKL_INTEGER64)
        SET(ARCH_PREFIX "_lp64")
    endif(MKL_INTEGER64)
endif(variable STREQUAL "32")


# Find include dir
find_path(MKL_INCLUDE_DIR mkl.h
    PATHS ${MKL_ROOT}/include)

# Find include directory
#  There is no include folder under linux
if(WIN32)
    find_path(INTEL_INCLUDE_DIR omp.h
        PATHS ${INTEL_ROOT}/include)
    set(MKL_INCLUDE_DIR ${MKL_INCLUDE_DIR} ${INTEL_INCLUDE_DIR})
endif()

# Find libraries

# Handle suffix
set(_MKL_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
#if(WIN32)
#    if(MKL_STATAIC)
#        set(CMAKE_FIND_LIBRARY_SUFFIXES .lib)
#    else()
#        set(CMAKE_FIND_LIBRARY_SUFFIXES _dll.lib)
#    endif()
#else()
#    if(MKL_STATAIC)
#        set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
#    else()
#        set(CMAKE_FIND_LIBRARY_SUFFIXES .so)
#    endif()
#endif()


# MKL is composed by four layers: Interface, Threading, Computational and RTL

if(MKL_SDL)
    find_library(MKL_LIBRARY mkl_rt
        PATHS ${MKL_ROOT}/lib/${MPL_ARCH}/)

    set(MKL_MINIMAL_LIBRARY ${MKL_LIBRARY})
else()
    ######################### Interface layer #######################
    if(WIN32)
        set(MKL_INTERFACE_LIBNAME mkl_intel_c)
    else()
        set(MKL_INTERFACE_LIBNAME mkl_intel_lp64) #ilp64
    endif()

    find_library(MKL_INTERFACE_LIBRARY ${MKL_INTERFACE_LIBNAME}
        PATHS ${MKL_ROOT}/lib/${MPL_ARCH}/)

    ######################## Threading layer ########################
    if(MKL_MULTI_THREADED)
        set(MKL_THREADING_LIBNAME mkl_intel_thread)
    else()
        set(MKL_THREADING_LIBNAME mkl_sequential)
    endif()

    find_library(MKL_THREADING_LIBRARY ${MKL_THREADING_LIBNAME}
        PATHS ${MKL_ROOT}/lib/${MPL_ARCH}/)

    ####################### Computational layer #####################
    find_library(MKL_CORE_LIBRARY mkl_core
        PATHS ${MKL_ROOT}/lib/${MPL_ARCH}/)
    find_library(MKL_FFT_LIBRARY mkl_cdft_core
        PATHS ${MKL_ROOT}/lib/${MPL_ARCH}/)
    find_library(MKL_SCALAPACK_LIBRARY mkl_scalapack_lp64     # _core / _ilp64
        PATHS ${MKL_ROOT}/lib/${MPL_ARCH}/)
    #find_library(MKL_BLACS_LIBRARY mkl_blacs_intelmpi_lp64     # / _ilp64, intelmpi, openmpi
    find_library(MKL_BLACS_LIBRARY mkl_blacs_lp64     # / _ilp64, intelmpi, openmpi
        PATHS ${MKL_ROOT}/lib/${MPL_ARCH}/)

    #find_library(MKL_THREAD_LIBRARY mkl_intel_thread
    #find_library(MKL_THREAD_LIBRARY mkl_pgi_thread
    #find_library(MKL_THREAD_LIBRARY mkl_gnu_thread
    #    PATHS ${MKL_ROOT}/lib/${MPL_ARCH}/)
    find_library(MKL_THREAD_LIBRARY pthread)

    find_library(MPI_LIBRARY mpi_mt
        PATHS /opt/intel/composer_xe_2013.0.079/mpirt/lib/${MPL_ARCH}/)

#message(blacs: "${MKL_ROOT}/lib/${MPL_ARCH}")
#set( MKL_BLACS_LIBRARY ${MKL_ROOT}/lib/${MPL_ARCH}/libmkl_blacs_lp64.a)

    ############################ RTL layer ##########################
    if(WIN32)
        set(MKL_RTL_LIBNAME iomp5md)
    else()
        set(MKL_RTL_LIBNAME iomp5)
    endif()
    find_library(MKL_RTL_LIBRARY ${MKL_RTL_LIBNAME}
        PATHS ${INTEL_ROOT}/lib/${MPL_ARCH} ${INTEL_MKL_ROOT}/lib )

set(MKL_LIBRARY ${MKL_INTERFACE_LIBRARY} 
		${MKL_THREADING_LIBRARY} 
		${MKL_CORE_LIBRARY} 
		${MKL_FFT_LIBRARY} 
		${MKL_SCALAPACK_LIBRARY} 
		${MKL_BLACS_LIBRARY} 
		${MKL_RTL_LIBRARY}
		${MPI_LIBRARY}
		${MKL_THREAD_LIBRARY}
		)

set(MKL_MINIMAL_LIBRARY ${MKL_INTERFACE_LIBRARY} 
			${MKL_THREADING_LIBRARY} 
			${MKL_CORE_LIBRARY} 
			${MKL_RTL_LIBRARY}
			)
endif()

#?
#set(CMAKE_FIND_LIBRARY_SUFFIXES ${_MKL_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})

find_package_handle_standard_args(MKL DEFAULT_MSG
    MKL_INCLUDE_DIR MKL_LIBRARY MKL_MINIMAL_LIBRARY)

if(MKL_FOUND)
    set(MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR})
    set(MKL_LIBRARIES ${MKL_LIBRARY})
    set(MKL_MINIMAL_LIBRARIES ${MKL_LIBRARY})
endif()

message("MKL_LIBRARIES: ${MKL_LIBRARIES}")
