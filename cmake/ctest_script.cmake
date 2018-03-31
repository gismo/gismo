######################################################################
## ctest_script.txt
## This file is part of the G+Smo library.
## https://raw.githubusercontent.com/gismo/gismo/stable/cmake/ctest_script.cmake
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012-2018
######################################################################

######################################################################
##
## To execute:
##
##   ctest -S ctest_script.cmake
##
## The script creates sources and build folders in the same directory
##
## For more verbosity add the flag
##
##   ctest -S /path/to/ctest_script.cmake -V
##
## or even -VV.
##
## Set execution options in the Configuration part.
## For multiple tests (eg. different compilers) make multiple copies
## of this file and adjust options. Few options can be passed by arguments:
##
## ctest -S ctest_script.cmake,"Experimental;Release;8;Ninja"
##
## On linux this script can be invoked in a cronjob. e.g.:
##    $ crontab -e
## Add the line:
##    0 3 * * * ctest -S /path/toctest_script.cmake,"Nightly" &>/dev/null
## save and exit. Now with
##    $ crontab -l
## you can see the scheduled task. The script will
## be executed every night at 03:00am.
##
######################################################################

## #################################################################
## Configuration
## #################################################################

# Test model (Nightly, Continuous, Experimental)
set(CTEST_TEST_MODEL Experimental)

# Configuration type (Debug Release RelWithDebInfo MinSizeRel)
set(CTEST_CONFIGURATION_TYPE Release)

# Number of jobs for build/test
#set(CTEST_BUILD_JOBS 8)
#set(CTEST_TEST_JOBS 10)

# The Generator for CMake
# ("Unix Makefiles", "Ninja", "Xcode", "NMake Makefiles", "NMake Makefiles JOM",
#  "MinGW Makefiles", "Visual Studio 12 2013", "Visual Studio 14 2015",
#  "Visual Studio 14 2015 Win64", and so on)
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")

set(CNAME cc)
set(CXXNAME g++)

# The above parameters can be reset by passing upto 6 arguments
# e.g. as: ctest -S ctest_script.cmake,"Experimental;Release;8;Ninja"
macro(read_args)
  set(narg ${ARGC})
    if (narg GREATER 0)
      set(CTEST_TEST_MODEL ${ARGV0})
    endif()
  if (narg GREATER 1)
    set(CTEST_CONFIGURATION_TYPE ${ARGV1})
  endif()
  if (narg GREATER 2)
    set(CTEST_BUILD_JOBS "${ARGV2}")
  endif()
  if (narg GREATER 3)
    set(CTEST_CMAKE_GENERATOR "${ARGV3}")
  endif()
  if (narg GREATER 5)
    set(CNAME "${ARGV4}")
    set(CXXNAME "${ARGV5}")
  endif()
endmacro(read_args)
read_args(${CTEST_SCRIPT_ARG})

# C/C++ compilers,  e.g. "cc/g++", "icc/icpc", "clang/clang++"
find_program (CC NAMES ${CNAME})
set(ENV{CC}  ${CC})
find_program (CXX NAMES ${CXXNAME})
set(ENV{CXX}  ${CXX})

# Other Environment variables and scripts
#set(ENV{CXXFLAGS} "-Ofast")
#execute_process(COMMAND source "/path/to/iccvars.sh intel64")
#set(ENV{LD_LIBRARY_PATH} /path/to/vendor/lib)
#set(ENV{MAKEFLAGS} "-j12")

# Build options
set(gismo_build_options
    -DGISMO_WARNINGS=OFF
    -DGISMO_COEFF_TYPE=double
    -DGISMO_BUILD_LIB=ON
    #-DCMAKE_CXX_STANDARD=11
    -DGISMO_BUILD_EXAMPLES=ON
    -DGISMO_BUILD_UNITTESTS=OFF
    #-DGISMO_WITH_OPENMP=ON
    #-DGISMO_WITH_MPI=ON
    #-DGISMO_WITH_SPECTRA=ON
    #-DGISMO_WITH_IPOPT=ON -DIpOpt_DIR=/path/to/ipopt
    #-DGISMO_WITH_PSOLID=ON -DParasolid_DIR=/path/to/parasolid
    #-DGISMO_BUILD_AXL=ON -DAxel_DIR=/path/to/axel
    -DGISMO_WITH_ONURBS=ON
    -DGISMO_WITH_TRILINOS=OFF
    -DGISMO_WITH_SPECTRA=OFF
    -DGISMO_EXTRA_DEBUG=OFF
    -DGISMO_BUILD_PCH=OFF
    #-DGISMO_PLAINDOX=ON
)
  
# Source folder (defaults inside the script directory)
set(CTEST_SOURCE_DIRECTORY ${CTEST_SCRIPT_DIRECTORY}/gismo_src)

# Build folder (defaults inside the script directory)
set(CTEST_BINARY_DIRECTORY ${CTEST_SCRIPT_DIRECTORY}/build_${CTEST_TEST_MODEL}_${CTEST_CONFIGURATION_TYPE}_${CXXNAME})

# Empty previous directory before building (otherwise builds are incremental)
#ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})
file(REMOVE_RECURSE ${CTEST_BINARY_DIRECTORY}/bin)

# Test timeout in seconds
set(CTEST_TEST_TIMEOUT 200)

# Coverage analysis
set(test_coverage FALSE)
if (test_coverage)
  find_program(CTEST_COVERAGE_COMMAND NAMES gcov)
  set(CTEST_CUSTOM_COVERAGE_EXCLUDE "${CTEST_SOURCE_DIRECTORY}/external/")
  set(ENV{CXXFLAGS} "$ENV{CXXFLAGS} -g -O0 --coverage -fprofile-arcs -ftest-coverage")
  set(ENV{CFLAGS} "$ENV{CFLAGS} -g -O0 --coverage -fprofile-arcs -ftest-coverage")
endif()

# Dynamic analysis
#Valgrind, Purify, BoundsChecker. ThreadSanitizer, AddressSanitizer,
#LeakSanitizer, MemorySanitizer, and UndefinedBehaviorSanitizer.
set(CTEST_MEMORYCHECK_TYPE "")

if ("${CTEST_MEMORYCHECK_TYPE}" STREQUAL "Valgrind")
  find_program(CTEST_MEMORYCHECK_COMMAND NAMES valgrind)
  set(MEMORYCHECK_SUPPRESSIONS_FILE "${CTEST_SOURCE_DIRECTORY}/cmake/valgrind_supp.txt")
  set(MEMORYCHECK_COMMAND_OPTIONS "--error-exitcode=1 --leak-check=yes -q")
  #--tool=memcheck --show-reachable=yes --num-callers=50 --track-origins=yes --trace-children=yes
elif(CTEST_MEMORYCHECK_TYPE STREQUAL "AddressSanitizer")
  set(ENV{CXXFLAGS} "$ENV{CXXFLAGS} -fsanitize=address -fno-omit-frame-pointer")
  set(ENV{LDFLAGS}  "$ENV{LDFLAGS} -fsanitize=address")
endif()

# Update type (eg. svn or git)
set(UPDATE_TYPE git)
find_program(CTEST_GIT_COMMAND NAMES ${UPDATE_TYPE})
set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")

# For continuous builds, number of seconds to stay alive
set(test_runtime 43200) #12h by default

# Computer ID shown on the dashboard (will be set automatically)
# set(CTEST_SITE "Server0407")

# Ignore certain tests during test or memcheck
#set(CTEST_CUSTOM_TESTS_IGNORE "")
#set(CTEST_CUSTOM_MEMCHECK_IGNORE "")

## #################################################################
## Test routines
## #################################################################

#message(STATUS "Preserve full output (CTEST_FULL_OUTPUT)")

# Initial checkout
if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
  #message("Initial checkout...")
  set(GISMO_REPOSITORY https://github.com/gismo/gismo.git)
  set(CTEST_CHECKOUT_COMMAND "${CTEST_GIT_COMMAND} clone --depth 1 --branch stable ${GISMO_REPOSITORY} gismo_src")
endif()

if("${CTEST_CMAKE_GENERATOR}" MATCHES "Make" OR "${CTEST_CMAKE_GENERATOR}" MATCHES "Ninja")
 set(CTEST_USE_LAUNCHERS 1)
else()
  set(CTEST_USE_LAUNCHERS 0)
endif()

set(ENV{CTEST_OUTPUT_ON_FAILURE} 1)
set( $ENV{LC_MESSAGES}      "en_EN" )

if(NOT DEFINED CTEST_TEST_MODEL AND DEFINED ENV{CTEST_TEST_MODEL})
  set(CTEST_TEST_MODEL $ENV{CTEST_TEST_MODEL})
endif()
if(NOT DEFINED CTEST_CONFIGURATION_TYPE AND DEFINED ENV{CTEST_CONFIGURATION_TYPE})
  set(CTEST_CONFIGURATION_TYPE $ENV{CTEST_CONFIGURATION_TYPE})
endif()

if(NOT DEFINED CTEST_SITE)
find_program(HOSTNAME_CMD NAMES hostname)
execute_process(COMMAND ${HOSTNAME_CMD} OUTPUT_VARIABLE HOSTNAME OUTPUT_STRIP_TRAILING_WHITESPACE)
set(CTEST_SITE "${HOSTNAME}")
endif()

# Name of this build
if(NOT DEFINED CTEST_BUILD_NAME)
find_program(UNAME NAMES uname)
execute_process(COMMAND "${UNAME}" "-s" OUTPUT_VARIABLE osname OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND "${UNAME}" "-m" OUTPUT_VARIABLE "cpu" OUTPUT_STRIP_TRAILING_WHITESPACE)
set(CTEST_BUILD_NAME "${osname}-${cpu} ${CTEST_CMAKE_GENERATOR}-${CTEST_CONFIGURATION_TYPE}-${CXXNAME}")
endif()

if(NOT CTEST_BUILD_JOBS)
include(ProcessorCount)
ProcessorCount(NPROC)
#message("Number of processors: ${NPROC}")
if(${NPROC} EQUAL 0)
  set(NPROC 1)
endif()
if(${NPROC} GREATER 20)
  set(CTEST_BUILD_JOBS 20)
else()
  math(EXPR CTEST_BUILD_JOBS "(1+${NPROC})>>1")
  #message("CTEST_BUILD_JOBS ${CTEST_BUILD_JOBS}")
endif()
endif()

if(NOT DEFINED CTEST_TEST_JOBS)
set(CTEST_TEST_JOBS ${CTEST_BUILD_JOBS})
endif()

if(${CTEST_CMAKE_GENERATOR} MATCHES "Unix Makefiles"
  OR "${CTEST_CMAKE_GENERATOR}" MATCHES "Ninja")
  set(CTEST_BUILD_FLAGS "-j ${CTEST_BUILD_JOBS}")
#message("Build flags: ${CTEST_BUILD_FLAGS}")
endif()

macro(run_ctests)
  ctest_configure(OPTIONS "${gismo_build_options};-DCTEST_USE_LAUNCHERS=${CTEST_USE_LAUNCHERS};-DDART_TESTING_TIMEOUT=${CTEST_TEST_TIMEOUT})")
  ctest_submit(PARTS Configure Update)
  ctest_build(TARGET gsUnitTest APPEND) # for older versions of ninja
  ctest_submit(PARTS Build)
  ctest_build(APPEND)
  ctest_submit(PARTS Build)
  ctest_build(TARGET unittests APPEND)
  ctest_submit(PARTS Build)
  ctest_test(PARALLEL_LEVEL ${CTEST_TEST_JOBS})
  ctest_submit(PARTS Test)

  if(test_coverage)
     #message("Running coverage..")
     ctest_coverage(BUILD "${CTEST_BINARY_DIRECTORY}" APPEND)
     ctest_submit(PARTS Coverage)
  endif()

  if(NOT "x${CTEST_MEMORYCHECK_TYPE}" STREQUAL "x")
    #message("Running memcheck..")
    ctest_memcheck()
    ctest_submit(PARTS MemCheck)
  endif()
endmacro(run_ctests)

file(MAKE_DIRECTORY "${CTEST_BINARY_DIRECTORY}")

ctest_start(${CTEST_TEST_MODEL})

if(NOT "${CTEST_TEST_MODEL}" STREQUAL "Continuous")

ctest_update()
run_ctests()

else() #continuous model

while(${CTEST_ELAPSED_TIME} LESS ${test_runtime})
  set(START_TIME ${CTEST_ELAPSED_TIME})
  ctest_update()
  if( ${count} GREATER 0 )
    run_ctests()
  endif()
  ctest_sleep(${START_TIME} 300 ${CTEST_ELAPSED_TIME})
endwhile()

endif(NOT "${CTEST_TEST_MODEL}" STREQUAL "Continuous")

# Cleanup xml files after upload
file(REMOVE_RECURSE ${CTEST_BINARY_DIRECTORY}/Testing)
