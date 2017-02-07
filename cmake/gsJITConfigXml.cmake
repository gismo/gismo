######################################################################
## gsJITConfigXml ---
## This file is part of the G+Smo library. 
##
## Author: Matthias Moeller
## Copyright (C) 2012 - 2017 RICAM-Linz.
######################################################################

## This macro defines the following variables:
##
## For <LANG> is C, CXX, CUDA, and Fortran
## JIT_<LANG>_COMPILER : JIT compiler command
## JIT_<LANG>_FLAGS    : JIT compiler flags
##
## JIT_INCLUDE_DIRS    : JIT compiler include directories

macro(gsJITConfigXml source_file target_file)

  # Gather information of programming languages not enabled by defualt
  # with G+Smo (CUDA and Fortran)
  include(CheckLanguage)
  check_language(CUDA)
  check_language(Fortran)
  
  # Set JIT compiler command
  set(JIT_C_COMPILER       ${CMAKE_C_COMPILER})  
  set(JIT_CXX_COMPILER     ${CMAKE_CXX_COMPILER})
  set(JIT_CUDA_COMPILER    ${CUDA_NVCC_EXECUTABLE})
  set(JIT_Fortran_COMPILER ${CMAKE_Fortran_COMPILER})

  # Set JIT compiler flags (build-type dependent)
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(JIT_C_FLAGS       ${CMAKE_C_FLAGS_DEBUG})
    set(JIT_CXX_FLAGS     ${CMAKE_CXX_FLAGS_DEBUG})
    set(JIT_CUDA_FLAGS    ${CUDA_NVCC_FLAGS_DEBUG})
    set(JIT_Fortran_FLAGS ${CMAKE_Fortran_FLAGS_DEBUG})
  elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(JIT_C_FLAGS       ${CMAKE_C_FLAGS_RELEASE})
    set(JIT_CXX_FLAGS     ${CMAKE_CXX_FLAGS_RELEASE})
    set(JIT_CUDA_FLAGS    ${CUDA_NVCC_FLAGS_RELEASE})
    set(JIT_Fortran_FLAGS ${CMAKE_Fortran_FLAGS_RELEASE})
  elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(JIT_C_FLAGS       ${CMAKE_C_FLAGS_RELWITHDEBINFO})
    set(JIT_CXX_FLAGS     ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
    set(JIT_CUDA_FLAGS    ${CUDA_NVCC_FLAGS_RELWITHDEBINFO})
    set(JIT_Fortran_FLAGS ${CMAKE_Fortran_FLAGS_RELWITHDEBINFO})
  elseif(CMAKE_BUILD_TYPE STREQUAL "MinSizeRel")
    set(JIT_C_FLAGS       ${CMAKE_C_FLAGS_MINSIZEREL})
    set(JIT_CXX_FLAGS     ${CMAKE_CXX_FLAGS_MINSIZEREL})
    set(JIT_CUDA_FLAGS    ${CUDA_NVCC_FLAGS_MINSIZEREL})
    set(JIT_Fortran_FLAGS ${CMAKE_Fortran_FLAGS_MINSIZEREL})
  endif()

  # Add JIT compiler flags (all build types)
  set(JIT_C_FLAGS       "${JIT_C_FLAGS} ${CMAKE_C_FLAGS} ${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS}")
  set(JIT_CXX_FLAGS     "${JIT_CXX_FLAGS} ${CMAKE_CXX_FLAGS} ${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS}")
  set(JIT_CUDA_FLAGS    "${JIT_CUDA_FLAGS} ${CUDA_NVCC_FLAGS} ${CMAKE_SHARED_LIBRARY_CREATE_CUDA_FLAGS}")
  set(JIT_Fortran_FLAGS "${JIT_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS} ${CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS}")

  # Add C++ standard
  if(CMAKE_CXX_STANDARD EQUAL "98")
    set(JIT_CXX_FLAGS "${JIT_CXX_FLAGS} ${CMAKE_CXX98_STANDARD_COMPILE_OPTION}")
    set(JIT_CUDA_FLAGS "${JIT_CXX_FLAGS} ${CUDA_CXX98_STANDARD_COMPILE_OPTION}")
  elseif(CMAKE_CXX_STANDARD EQUAL "11")
    set(JIT_CXX_FLAGS "${JIT_CXX_FLAGS} ${CMAKE_CXX11_STANDARD_COMPILE_OPTION}")
    set(JIT_CUDA_FLAGS "${JIT_CXX_FLAGS} ${CUDA_CXX11_STANDARD_COMPILE_OPTION}")
  elseif(CMAKE_CXX_STANDARD EQUAL "14")
    set(JIT_CXX_FLAGS "${JIT_CXX_FLAGS} ${CMAKE_CXX14_STANDARD_COMPILE_OPTION}")
    set(JIT_CUDA_FLAGS "${JIT_CXX_FLAGS} ${CUDA_CXX14_STANDARD_COMPILE_OPTION}")
  endif()
  
  # Set includes
  string(REPLACE ";" " -I" JIT_INCLUDES "-I${GISMO_DEV_INCLUDE_DIRS} -I${GISMO_INCLUDE_DIRS} -I${CMAKE_CURRENT_SOURCE_DIR} -I${CMAKE_CURRENT_SOURCE_DIR}/src")

  # Set library directories
  string(REPLACE ";" " -L" JIT_LIBRARIES "-L${GISMO_LIBRARY_DIR} -L${GISMO_DEV_LIBRARY_DIR}")

  # Set libraries
  string(REPLACE ";" " -l" JIT_LIBRARIES "${JIT_LIBRARIES} -l${GISMO_LIBRARIES} -l${GISMO_DEV_LIBRARIES}")
  
  # Set extra libraries
  if (DEFINED gismo_LINKER)
    string(REPLACE ";" " -l" JIT_LIBRARIES "${JIT_LIBRARIES} -l${gismo_LINKER}")
  endif()
  
  # Generate config XML file
  configure_file(${source_file} ${target_file})
  
endmacro()