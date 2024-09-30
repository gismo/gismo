######################################################################
## gsJITConfigXml.cmake
## This file is part of the G+Smo library. 
##
## Author: Matthias Moeller
######################################################################

## This macro defines the following variables:
##
## For <LANG> is C, CXX, CUDA, and Fortran
## JIT_<LANG>_COMPILER : JIT compiler command
## JIT_<LANG>_FLAGS    : JIT compiler flags
##
## JIT_INCLUDE_DIRS    : JIT compiler include directories

macro(gsJITConfigXml source_file target_file)

  set(extra_macro_args ${ARGN})
  list(LENGTH extra_macro_args num_extra_macro_args)

  if(${num_extra_macro_args} GREATER 0)
    
    # Gather information of programming languages not enabled by
    # default with G+Smo (CUDA and Fortran)
    include(CheckLanguage)

    foreach(extra_macro_arg ${extra_macro_args})
      if("${extra_macro_arg}" STREQUAL "CUDA")
        # Check for CUDA 
        check_language(CUDA)
        if(CMAKE_CUDA_COMPILER)
          enable_language(Fortran)
        endif()
      elseif("${extra_macro_arg}" STREQUAL "FORTRAN")
        # Check for Fortran
        check_language(Fortran)
        if(CMAKE_Fortran_COMPILER)
          enable_language(Fortran)
        endif()
      else()
        message(FATAL_ERROR "Unsupport optional argument ${extra_macro_arg}.")
      endif()
    endforeach()
  endif()
  
  # Set JIT compiler command
  set(JIT_C_COMPILER       ${CMAKE_C_COMPILER})  
  set(JIT_CXX_COMPILER     ${CMAKE_CXX_COMPILER})
  set(JIT_CUDA_COMPILER    ${CUDA_NVCC_EXECUTABLE})
  set(JIT_Fortran_COMPILER ${CMAKE_Fortran_COMPILER})

  # Get build-type as upper-case string
  string(TOUPPER ${CMAKE_BUILD_TYPE} JIT_BUILD_TYPE)
  
  # Set JIT compiler flags (build-type dependent)
  set(JIT_C_FLAGS       ${CMAKE_C_FLAGS_${JIT_BUILD_TYPE}})
  set(JIT_CXX_FLAGS     ${CMAKE_CXX_FLAGS_${JIT_BUILD_TYPE}})
  set(JIT_CUDA_FLAGS    ${CUDA_NVCC_FLAGS_${JIT_BUILD_TYPE}})
  set(JIT_Fortran_FLAGS ${CMAKE_Fortran_FLAGS_${JIT_BUILD_TYPE}})

  # Add JIT compiler flags (all build types)
  set(JIT_C_FLAGS "${JIT_C_FLAGS} ${CMAKE_C_FLAGS} ${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS}")
  string(REPLACE "-fvisibility=hidden"         "" JIT_C_FLAGS ${JIT_C_FLAGS})
  string(REPLACE "-fvisibility-inlines-hidden" "" JIT_C_FLAGS ${JIT_C_FLAGS})

  # Create a set of shared library variable specific to C++
  # For 90% of the systems, these are the same flags as the C versions
  # so if these are not set just copy the flags from the c version
  if(NOT DEFINED CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS)
    set(JIT_CXX_FLAGS "${JIT_CXX_FLAGS} ${CMAKE_CXX_FLAGS} ${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS}")
  else()
    set(JIT_CXX_FLAGS "${JIT_CXX_FLAGS} ${CMAKE_CXX_FLAGS} ${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS}")
  endif()

  # Create a set of shared library variable specific to CUDA
  # For 90% of the systems, these are the same flags as the C versions
  # so if these are not set just copy the flags from the c version
  if(NOT DEFINED CMAKE_SHARED_LIBRARY_CREATE_CUDA_FLAGS)
  else()
    set(JIT_CUDA_FLAGS "${JIT_CUDA_FLAGS} ${CUDA_NVCC_FLAGS} ${CMAKE_SHARED_LIBRARY_CREATE_CUDA_FLAGS}")
  endif()
  
  # Create a set of shared library variable specific to Fortran
  # For 90% of the systems, these are the same flags as the C versions
  # so if these are not set just copy the flags from the c version
  if(NOT DEFINED CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS)
    set(JIT_Fortran_FLAGS "${JIT_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS} ${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS}")
  else()
    set(JIT_Fortran_FLAGS "${JIT_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS} ${CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS}")
  endif()
  
  # Add C++ standard and pic
  set(JIT_CXX_FLAGS "${JIT_CXX_FLAGS} ${CMAKE_CXX${CMAKE_CXX_STANDARD}_STANDARD_COMPILE_OPTION} ${CMAKE_CXX_COMPILE_OPTIONS_PIC}")
  set(JIT_CUDA_FLAGS "${JIT_CXX_FLAGS} ${CUDA_CXX${CMAKE_CXX_STANDARD}_STANDARD_COMPILE_OPTION}")
  
  # PIC for 64bit machines
  #if (CMAKE_COMPILER_IS_GNUCXX AND CMAKE_SIZEOF_VOID_P EQUAL 8)
  #    set(JIT_CXX_FLAGS "${JIT_CXX_FLAGS} ${CMAKE_CXX_COMPILE_OPTIONS_PIC}")
  #endif()

  string(REPLACE "-fvisibility=hidden"         "" JIT_CXX_FLAGS ${JIT_CXX_FLAGS})
  string(REPLACE "-fvisibility-inlines-hidden" "" JIT_CXX_FLAGS ${JIT_CXX_FLAGS})

  # Set includes
  string(REPLACE ";" " -I" JIT_INCLUDES "-I${GISMO_INCLUDE_DIRS} -I${CMAKE_CURRENT_SOURCE_DIR} -I${CMAKE_CURRENT_SOURCE_DIR}/src  -I${gismo_BINARY_DIR}")

  # Set (extra) library directories
  if(GISMO_LIBRARY_DIR)
    string(REPLACE ";" " -L" JIT_LIBRARIES "-L${GISMO_LIBRARY_DIR}")
  endif()

  foreach(lib ${GISMO_LIBRARIES} ${gismo_LINKER})
    if (NOT IS_ABSOLUTE ${lib})
      set(JIT_LIBRARIES "${JIT_LIBRARIES} -l${lib}")
    else()
      set(JIT_LIBRARIES "${JIT_LIBRARIES} ${lib}")
    endif()
  endforeach()
  
  # Set libraries
  #string(REPLACE ";" " -l" JIT_LIBRARIES "${JIT_LIBRARIES} -l${GISMO_LIBRARIES}")
  
  # Set extra libraries
  #if (DEFINED gismo_LINKER)
  #  string(REPLACE ";" " -l" JIT_LIBRARIES "${JIT_LIBRARIES} -l${gismo_LINKER}")
  #endif()
  
  # Generate config XML file
  configure_file(${source_file} ${target_file})
  
endmacro()
