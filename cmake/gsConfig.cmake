######################################################################
## CMakeLists.txt ---
## This file is part of the G+Smo library. 
##
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - 2015 RICAM-Linz.
######################################################################

## #################################################################
## Configuration
## #################################################################

# Set a default coefficient numeric types if not specified
if(NOT GISMO_COEFF_TYPE)
  set (GISMO_COEFF_TYPE "double" CACHE STRING
   "Coefficient type(float, double, long double, mpfr::mpreal)" FORCE)
   set_property(CACHE GISMO_COEFF_TYPE PROPERTY STRINGS
   "float" "double" "long double" "mpfr::mpreal"
   )
endif()

if(NOT GISMO_INDEX_TYPE)
  set (GISMO_INDEX_TYPE "int" CACHE STRING
   "Index type(int, unsigned, size_t)" FORCE)
   set_property(CACHE GISMO_INDEX_TYPE PROPERTY STRINGS
   "int" "unsigned" "size_t"
   )
endif()

# Shared pointer
find_package (TR1 QUIET)

## #################################################################
## Setup build types
## #################################################################

SET( CMAKE_CXX_FLAGS_MAINTAINER "-Wall -Wabi -DUSE_GISMO_STACK_WALKER"
     CACHE STRING
    "Flags used by the C++ compiler during maintainer builds."
    FORCE )
SET( CMAKE_C_FLAGS_MAINTAINER "-Wall -pedantic" CACHE STRING
    "Flags used by the C compiler during maintainer builds."
    FORCE )
SET( CMAKE_EXE_LINKER_FLAGS_MAINTAINER
    "-Wl,--warn-unresolved-symbols,--warn-once" CACHE STRING
    "Flags used for linking binaries during maintainer builds."
    FORCE )
SET( CMAKE_SHARED_LINKER_FLAGS_MAINTAINER
    "-Wl,--warn-unresolved-symbols,--warn-once" CACHE STRING
    "Flags used by the shared libraries linker during maintainer builds."
    FORCE )
MARK_AS_ADVANCED(
    CMAKE_CXX_FLAGS_MAINTAINER
    CMAKE_C_FLAGS_MAINTAINER
    CMAKE_EXE_LINKER_FLAGS_MAINTAINER
    CMAKE_SHARED_LINKER_FLAGS_MAINTAINER )

# Update the documentation string of CMAKE_BUILD_TYPE for GUIs
SET( CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel Maintainer."
    FORCE )

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
   #set(CMAKE_BUILD_TYPE Debug CACHE STRING 
   set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING 
   "Type of build (Debug, Release, RelWithDebInfo, MinSizeRel)" FORCE)
   # Set the possible values of build type for cmake-gui
   set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
     "MinSizeRel" "RelWithDebInfo")
endif()


# Remove NDEBUG flag from RelWithDebInfo builds
if(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    STRING(REPLACE "-DNDEBUG" "" CMAKE_CXX_FLAGS_RELWITHDEBINFO ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
endif()

set(${PROJECT_NAME}_ARCHIVE_OUTPUT_DIRECTORY lib)
set(${PROJECT_NAME}_RUNTIME_OUTPUT_DIRECTORY bin)
set(${PROJECT_NAME}_LIBRARY_OUTPUT_DIRECTORY lib)
foreach(config ${CMAKE_CONFIGURATION_TYPES}) # For Visual studio
    # overrides Debug/Release subfolders
    string(TOUPPER ${config} CONFIG) 
    set(${PROJECT_NAME}_ARCHIVE_OUTPUT_DIRECTORY_${CONFIG} lib)
    set(${PROJECT_NAME}_RUNTIME_OUTPUT_DIRECTORY_${CONFIG} bin)
    set(${PROJECT_NAME}_LIBRARY_OUTPUT_DIRECTORY_${CONFIG} lib)
endforeach() 

#Enable C++ 11
if(GISMO_BUILD_CPP11 AND UNIX)
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
endif()

# Print compilation statistics (these flags work on GCC compiler only)
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftime-report")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Q")

if (GISMO_BUILD_COVERAGE AND CMAKE_COMPILER_IS_GNUCXX)
  # see http://www.cmake.org/Wiki/CTest:Coverage
  # and http://cmake.3232098.n2.nabble.com/Running-coverage-analysis-td7145452.html
  include(CodeCoverage)
  SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ftest-coverage -fprofile-arcs")
  SET(CMAKE_EXE_LINKER_FLAGS "-fprofile-arcs -ftest-coverage")
endif(GISMO_BUILD_COVERAGE AND CMAKE_COMPILER_IS_GNUCXX)

if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
include( OptimizeForArchitecture )
endif("${CMAKE_BUILD_TYPE}" STREQUAL "Release")

if("x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xMSVC")

    #  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 18)
    # Update 4 of MSVC 2013 ?
    #    endif()

    # Disable checked iterators
    set(CMAKE_CXX_FLAGS    "${CMAKE_CXX_FLAGS}  /bigobj /D_SECURE_SCL=0")
    # See http://msdn.microsoft.com/en-us/library/hh697468.aspx
    #add_definitions(-D_HAS_ITERATOR_DEBUGGING=0)
    #add_definitions(-D_SECURE_SCL=0)
    #add_definitions(-D_ITERATOR_DEBUG_LEVEL=0) #VS2012

    # disable incremental linking for executables (it doesn't help for linking with libraries) -- check
    STRING(REPLACE "/INCREMENTAL:YES" "/INCREMENTAL:NO" CMAKE_EXE_LINKER_FLAGS_DEBUG ${CMAKE_EXE_LINKER_FLAGS_DEBUG})
    STRING(REPLACE "/INCREMENTAL:YES" "/INCREMENTAL:NO" CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO ${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO})

#    if ( GISMO_BUILD_LIB )
#    # /MD /MDd
#      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MD")
#    endif()

    if (CMAKE_SIZEOF_VOID_P EQUAL 8) #64bit compiler 
       # Note: On 64bit-platforms, /Wp64 flag is present, causing extra warnings
       set(CMAKE_CXX_FLAGS    "${CMAKE_CXX_FLAGS} /wd4244 /wd4267")
    #else() #32bit compiler has CMAKE_SIZEOF_VOID_P EQUAL 4
    endif()

endif()

if(GISMO_EXTRA_DEBUG)
  include(gsDebugExtra)
endif(GISMO_EXTRA_DEBUG)

if("x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xMSVC")
  # Force to always compile with W4
  if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
  else()
    set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS} /W4")
  endif()

elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # Update if necessary
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long") # -Woverloaded-virtual -Wconversion -Wextra -pedantic
  if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftrack-macro-expansion=0")
  endif()
endif()

if (MINGW)
  # export explicit template instantiations
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--export-all-symbols")

  # large files can overflow pe/coff sections, so use the pe+ format
  include(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-Wa,-mbig-obj" HAS_MBIGOBJ)
  if(NOT HAS_MBIGOBJ)
    message(WARNING "Current compiler does not suppport -Wa,-mbig-obj option.")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffunction-sections -Wl,--gc-sections")
  endif()
endif()

if (GISMO_WITH_OPENMP)
   find_package(OpenMP REQUIRED)
   set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if (CMAKE_COMPILER_IS_GNUCXX AND NOT GISMO_WITH_OPENMP)
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
endif()

#message("CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS}")
#message("CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG}")
#message("CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE}")
#message("CMAKE_CXX_FLAGS_RELWITHDEBINFO ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
#string(TOUPPER ${CMAKE_BUILD_TYPE} TEMP)
#message(STATUS "Using compilation flags: ${CMAKE_CXX_FLAGS}, ${CMAKE_CXX_FLAGS_${TEMP}}")
