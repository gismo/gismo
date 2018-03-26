######################################################################
## ctest_script.cmake
## This file is part of the G+Smo library. 
##
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - 2018 RICAM-Linz.
##
## To execute:
##
##   ctest -j jobs -S /path/to/ctest_script.cmake
##
## It is recommended to make a copy of the file (especially using git).
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
## of this file and adjust options.
##
## On linux this script can be invoked in a cronjob. e.g.:
##    $ >crontab -e
## Add the line:
##    0 3 * * * /path/to/script/nightly_cron.sh &>/dev/null
## save and exit. Now with
##    $ crontab -l
## you can see the scheduled task.
##
## "0 3 * * * " means that the script will be executed
## every night at 03:00am.
##
######################################################################

## #################################################################
## Configuration
## #################################################################

# Test type (Nightly, Continuous, Experimental)
set(test_model "Experimental")

# Build type
set(CTEST_CONFIGURATION_TYPE RelWithDebInfo)

# The Generator for CMake.
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
#set(CTEST_CMAKE_GENERATOR "Ninja")
#set(CTEST_CMAKE_GENERATOR "NMake Makefiles")
#set(CTEST_CMAKE_GENERATOR "NMake Makefiles JOM")
#set(CTEST_CMAKE_GENERATOR "MinGW Makefiles")
#set(CTEST_CMAKE_GENERATOR "Visual Studio 12 2013")
#set(CTEST_CMAKE_GENERATOR "Visual Studio 14 2015")
#set(CTEST_CMAKE_GENERATOR "Visual Studio 14 2015 Win64")
#set(CTEST_CMAKE_GENERATOR "Xcode")
#set(CTEST_CMAKE_GENERATOR "CodeBlocks")
#set(CTEST_CMAKE_GENERATOR "Sublime Text 2")
#set(CTEST_CMAKE_GENERATOR "Eclipse CDT4")
#set(CTEST_BUILD_FLAGS "-j12")

# Set environment/compiler
#set(ENV{MAKEFLAGS} "-j12")
#set(ENV{LD_LIBRARY_PATH} /path/to/vendor/lib)
#set(ENV{CC}  "gcc")
#set(ENV{CXX} "g++")
#exec_program(source ARGS "/path/to/iccvars.sh intel64")
#set(ENV{CC}  "icc")
#set(ENV{CXX} "icpc")
#set(ENV{CC}  "clang")
#set(ENV{CXX} "clang++")

# Build options
set(gismo_build_options 
    -DGISMO_BUILD_LIB=ON
    -DGISMO_WITH_ONURBS=ON
    #-DGISMO_WITH_IPOPT=ON
    #-DIpOpt_DIR=/path/to/ipopt
    #-DGISMO_WITH_PSOLID=ON
    #-DParasolid_DIR=/path/to/parasolid
    #-DGISMO_BUILD_AXL=ON -DAxel_DIR=/path/to/axel
    #-DGISMO_PLAINDOX=ON
    #-DGISMO_BUILD_COVERAGE=ON
    -DGISMO_BUILD_UNITTESTS=ON
)

# For continuous builds, number of seconds to stay alive
set(test_runtime 40000)

# ID for this computer that shows up on the dashboard.
if(NOT HOSTNAME)
find_program(HOSTNAME_CMD NAMES hostname)
exec_program(${HOSTNAME_CMD} ARGS OUTPUT_VARIABLE HOSTNAME)
set(CTEST_SITE "${HOSTNAME}")
endif(NOT HOSTNAME)

# Name of this build
if(NOT CTEST_BUILD_NAME)
find_program(UNAME NAMES uname)
exec_program("${UNAME}" ARGS "-s" OUTPUT_VARIABLE osname)
exec_program("${UNAME}" ARGS "-m" OUTPUT_VARIABLE "cpu")
set(CTEST_BUILD_NAME "${osname}-${cpu} ${CTEST_CMAKE_GENERATOR}/${CTEST_CONFIGURATION_TYPE} $ENV{CXX}")
endif(NOT CTEST_BUILD_NAME)
  
# Source folder (defaults inside the script directory)
set(CTEST_SOURCE_DIRECTORY ${CTEST_SCRIPT_DIRECTORY}/gismo_src)

#Build folder (defaults inside the script directory)
set(CTEST_BINARY_DIRECTORY ${CTEST_SCRIPT_DIRECTORY}/build_ctest)

# Update type (eg. svn or git)
set(UPDATE_TYPE git)
set(CTEST_UPDATE_COMMAND "git")
set(CTEST_GIT_COMMAND "git")

# Timeouts
set(CTEST_TEST_TIMEOUT 800 CACHE STRING 
    "Maximum time allowed before CTest will kill the test.") 
set(DART_TESTING_TIMEOUT 800 CACHE STRING 
    "Maximum time allowed before CTest will kill the test." FORCE)

# Coverage analysis
#set(test_coverage TRUE)
#set(CTEST_COVERAGE_COMMAND "/usr/bin/gcov")
#set(CTEST_CUSTOM_COVERAGE_EXCLUDE "${CTEST_SOURCE_DIRECTORY}/external/")

# Memory check with valgrind
#set(test_memcheck true)
#set(CTEST_MEMORYCHECK_COMMAND "/usr/bin/valgrind")
#set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE "/path_to/suppression_file.supp")
#set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --leak-check=full --show-reachable=yes --track-origins=yes")

# Initial checkout
set(GISMO_REPOSITORY https://github.com/gismo/gismo.git)
if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
  message("Initial checkout...")
  set(CTEST_CHECKOUT_COMMAND "git clone --depth 1 --branch stable ${GISMO_REPOSITORY} gismo_src")
endif()

## #################################################################
## Test routines
## #################################################################

macro(run_ctests)
  ctest_configure(OPTIONS "${gismo_build_options}")
  ctest_submit(PARTS Update Notes Configure)
  ctest_build(TARGET gsUnitTest) # for older versions of ninja
  ctest_submit(PARTS Build)
  ctest_build()
  ctest_build(TARGET unittests APPEND)
  ctest_submit(PARTS Build)
  ctest_test()
  ctest_submit(PARTS Test)

  if(test_coverage)
     message("Running coverage..")
     ctest_coverage(BUILD "${CTEST_BINARY_DIRECTORY}" APPEND)
     ctest_submit(PARTS Coverage)
  endif()

  if(test_memcheck)
    message("Running memcheck..")
    ctest_memcheck()
    ctest_submit(PARTS MemCheck)
  endif()
endmacro(run_ctests)

file(MAKE_DIRECTORY "${CTEST_BINARY_DIRECTORY}")

# Empty existing directory before building
#ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})

ctest_start(${test_model})

if(NOT "${test_model}" STREQUAL "Continuous")

ctest_update()
run_ctests()

else() #continuous model

while(${CTEST_ELAPSED_TIME} LESS ${test_runtime})
  set(START_TIME ${CTEST_ELAPSED_TIME})
  ctest_update(RETURN_VALUE count)
  #message(STATUS "Found ${count} changed files.")
  if( ${count} GREATER 0 )
    run_ctests()
  endif()
  ctest_sleep(${START_TIME} 300 ${CTEST_ELAPSED_TIME})
endwhile()

endif(NOT "${test_model}" STREQUAL "Continuous")
