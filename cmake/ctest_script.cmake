######################################################################
## ctest_script.cmake
## This file is part of the G+Smo library. 
##
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - 2016 RICAM-Linz.
######################################################################

## #################################################################
## Configuration
## #################################################################

# ID for this computer that shows up on the dashboard.
set(CTEST_SITE "SP2_0407_OpenSuse")

# Build type
set(CTEST_BUILD_CONFIGURATION RelWithDebInfo)

# Test type (Nightly, Continuous, Experimental)
set(dashboard_model "Nightly")

# Build flags
set( CTEST_BUILD_FLAGS "-j2")

# Source folder
set(CTEST_SOURCE_DIRECTORY ${CTEST_DASHBOARD_ROOT}/..)

#Build folder
set(CTEST_BINARY_DIRECTORY ${CTEST_SOURCE_DIRECTORY}/build_ctest)

# Name of this build (eg. compiler/build type)
set(CTEST_BUILD_NAME "GCC_4.7_${CTEST_BUILD_CONFIGURATION}")

# Update type (eg. svn or git)
set( UPDATE_TYPE svn)
#set( CTEST_UPDATE_COMMAND "svn_github.sh")

# Timeouts
set(CTEST_TEST_TIMEOUT 200 CACHE STRING 
    "Maximum time allowed before CTest will kill the test.") 
set(DART_TESTING_TIMEOUT 200 CACHE STRING 
    "Maximum time allowed before CTest will kill the test." FORCE)

# Build options
set(gismo_build_options 
    #-DGISMO_BUILD_COVERAGE=ON
    -DGISMO_PLAINDOX=ON # plain doxygen output for the trac Wiki
    -DGISMO_BUILD_LIB=ON
    -DGISMO_WITH_ONURBS=ON
    -DGISMO_WITH_IPOPT=ON
    -DIpOpt_DIR=/home/amantzaflaris/dashboards/IpOpt_inst
    -DGISMO_WITH_PSOLID=ON
    -DParasolid_DIR=/home/amantzaflaris/dropboxOUT/gforge/parasolid/26/1/153
    #-DGISMO_BUILD_AXL=ON -DAxel_DIR=/home/amantzaflaris/build/axel
)

# Coverage analysis
#set(dashboard_do_coverage TRUE)
#set(CTEST_COVERAGE_COMMAND "/usr/bin/gcov")
#set(CTEST_CUSTOM_COVERAGE_EXCLUDE "${CTEST_SOURCE_DIRECTORY}/external/")

# Memory check with valgrind
set(dashboard_do_memcheck true)
set(CTEST_MEMORYCHECK_COMMAND "/usr/bin/valgrind")

# The Generator for CMake.
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")

# This is the directory where the source and build trees will be placed.
get_filename_component(CTEST_DASHBOARD_ROOT "${CTEST_SCRIPT_DIRECTORY}/dashboard_gismo" ABSOLUTE)


## #################################################################
## Test routine
## #################################################################

if("${dashboard_model}" STREQUAL "Continuous")
  set(dashboard_continuous 1)
  set(dashboard_loop 40000)
else()
  set(dashboard_continuous 0)
  set(dashboard_loop 0)
endif()


#todo
#if("${dashboard_model}" STREQUAL "Continuous")

#if(NOT EXISTS "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt")
#  message("Starting fresh configuration...")
#  write_cache()
#endif()

file(MAKE_DIRECTORY "${CTEST_BINARY_DIRECTORY}")
ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})

ctest_start(${dashboard_model})
ctest_configure(OPTIONS "${gismo_build_options}")
ctest_submit(PARTS Update Notes Configure)
ctest_build()
ctest_submit(PARTS Build)
ctest_test()
ctest_submit(PARTS Test)

if(dashboard_do_coverage)
   message("Running coverage..")
   ctest_coverage(BUILD "${CTEST_BINARY_DIRECTORY}" APPEND)
   ctest_submit(PARTS Coverage)
endif()

if(dashboard_do_memcheck)
    message("Running memcheck..")
    ctest_memcheck()
    ctest_submit(PARTS MemCheck)
endif()

