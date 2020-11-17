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
## Options can be passed by arguments (default options are displayed here):
## NOTE: Line wrap can be done under Linux/macOS with \, under Windows with ^
##
## ctest -S ctest_script.cmake -D CTEST_TEST_MODEL=Experimental \
##   -D CTEST_CONFIGURATION_TYPE=Release -D CTEST_BUILD_JOBS=8 \
##   -D CTEST_CMAKE_GENERATOR="Unix Makefiles" -D CNAME=gcc -D CXXNAME=g++ \
##   -D CTEST_TEST_TIMEOUT=100 -D CTEST_MEMORYCHECK_TYPE=Valgrind \
##   -D DO_COVERAGE=TRUE
##
## Different dashboard projects and subprojects are possible:
##
## ctest -S /path/to/ctest_script.cmake -D PROJECT_NAME=myGismo \
##   -D CTEST_BUILD_JOBS=2 -D CTEST_CMAKE_GENERATOR="Unix Makefiles" \
##   -D CTEST_TEST_TIMEOUT=100 -D LABELS_FOR_SUBPROJECTS='gismo;examples;unittests' \
##   -D CTEST_SOURCE_DIRECTORY=./gismo_src -D CTEST_BINARY_DIRECTORY=./build
##
## On linux this script can be invoked in a cronjob. e.g.:
##    $ crontab -e
## Add the line:
##    0 3 * * * ctest -S /path/toctest_script.cmake -D CTEST_TEST_MODEL=Nightly -Q
## save and exit. Now with
##    $ crontab -l
## you can see the scheduled task. The script will
## be executed every night at 03:00am.
##
######################################################################

######################################################################
##
## Complete Options List:
##
## ctest parameters: (ctest -D ...)
##   CMAKE_ARGS
##   CNAME
##   CTEST_BINARY_DIRECTORY
##   CTEST_BUILD_JOBS
##   CTEST_BUILD_NAME
##   CTEST_CMAKE_GENERATOR
##   CTEST_CONFIGURATION_TYPE
##   CTEST_COVERAGE_COMMAND
##   CTEST_MEMORYCHECK_TYPE
##   CTEST_SITE
##   CTEST_SOURCE_DIRECTORY
##   CTEST_TEST_JOBS
##   CTEST_TEST_MODEL
##   CTEST_TEST_TIMEOUT
##   CXXNAME
##   DO_COVERAGE
##   DROP_LOCATION
##   DROP_METHOD
##   DROP_SITE
##   EMPTY_BINARY_DIRECTORY
##   GISMO_BRANCH
##   GISMO_SUBMODULES
##   LABELS_FOR_SUBPROJECTS
##   PROJECT_NAME
##   UPDATE_MODULES
##   UPDATE_TYPE
##
## Environment
##   CFLAGS
##   CXXFLAGS
##   LDFLAGS
##
######################################################################

cmake_minimum_required(VERSION 2.8.8)

if (POLICY CMP0048)# CMake 3.0
  cmake_policy(SET CMP0011 NEW)
  cmake_policy(SET CMP0042 NEW)
  cmake_policy(SET CMP0048 NEW)
endif()

if (POLICY CMP0054)# CMake 3.1
  cmake_policy(SET CMP0054 NEW)
endif()

if (POLICY CMP0053)# CMake 3.1.3
  cmake_policy(SET CMP0053 NEW)
endif()

if (POLICY CMP0063)# CMake 3.3
  cmake_policy(SET CMP0063 NEW)
endif()

## #################################################################
## Configuration
## #################################################################

# Test model (Nightly, Continuous, Experimental)
if (NOT DEFINED CTEST_TEST_MODEL)
  set(CTEST_TEST_MODEL Experimental)
endif()

# Configuration type (Debug Release RelWithDebInfo MinSizeRel)
if (NOT DEFINED CTEST_CONFIGURATION_TYPE)
  set(CTEST_CONFIGURATION_TYPE Release)
endif()

# Number of jobs for build/test (later on)
#set(CTEST_BUILD_JOBS 8)
#set(CTEST_TEST_JOBS 10)

# Tip for C/C++ compilers
# e.g. "cc/g++", "icc/icpc", "clang/clang++", "mpicc/mpic++", cl.exe/cl.exe
#set(CNAME cc)
#set(CXXNAME g++)

# Test timeout in seconds
if (NOT DEFINED CTEST_TEST_TIMEOUT)
  set(CTEST_TEST_TIMEOUT 200)
endif()

# Dynamic analysis
#Valgrind, Purify, BoundsChecker. ThreadSanitizer, AddressSanitizer,
#LeakSanitizer, MemorySanitizer, and UndefinedBehaviorSanitizer.
if (NOT DEFINED CTEST_MEMORYCHECK_TYPE)
  set(CTEST_MEMORYCHECK_TYPE "None")
endif()

# Coverage analysis - only GCC
# if GCC was changed with CNAME/CXXNAME, CTEST_COVERAGE_COMMAND needs
# also to be changed
if (NOT DEFINED DO_COVERAGE)
  set(DO_COVERAGE FALSE)
endif()

# The above parameters can be reset by passing upto 9 arguments
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
  if (narg GREATER 4)
    set(CNAME "${ARGV4}")
  endif()
  if (narg GREATER 5)
    set(CXXNAME "${ARGV5}")
  endif()
  if (narg GREATER 6)
    set(CTEST_TEST_TIMEOUT "${ARGV6}")
  endif()
  if (narg GREATER 7)
    set(CTEST_MEMORYCHECK_TYPE "${ARGV7}")
  endif()
  if (narg GREATER 8)
    set(DO_COVERAGE "${ARGV8}")
  endif()
endmacro(read_args)
read_args(${CTEST_SCRIPT_ARG})

if(DEFINED CNAME)
  find_program (CC NAMES ${CNAME})
  set(ENV{CC}  ${CC})
endif()
if(DEFINED CXXNAME)
  find_program (CXX NAMES ${CXXNAME})
  set(ENV{CXX}  ${CXX})
endif()

# Other Environment variables and scripts
#set(ENV{OMP_NUM_THREADS} 3)
#set(ENV{CXXFLAGS} "-Ofast")
#execute_process(COMMAND source "/path/to/iccvars.sh intel64")
#set(ENV{LD_LIBRARY_PATH} /path/to/vendor/lib)
#set(ENV{MAKEFLAGS} "-j12")

# Build options
if(NOT DEFINED CMAKE_ARGS)
  set(CMAKE_ARGS
    -DGISMO_WARNINGS=OFF
    -DGISMO_COEFF_TYPE=double
    -DGISMO_BUILD_LIB=ON
    #-DCMAKE_CXX_STANDARD=11
    -DGISMO_BUILD_EXAMPLES=ON
    -DGISMO_BUILD_UNITTESTS=ON
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
    -DNOSNIPPETS=OFF
    )
endif()

# Source folder (defaults inside the script directory)
if(NOT DEFINED CTEST_SOURCE_DIRECTORY)
  set(CTEST_SOURCE_DIRECTORY ${CTEST_SCRIPT_DIRECTORY}/gismo_src)
endif()

# Build folder (defaults inside the script directory)
if(NOT DEFINED CTEST_BINARY_DIRECTORY)
  set(CTEST_BINARY_DIRECTORY ${CTEST_SCRIPT_DIRECTORY}/build_${CTEST_TEST_MODEL}${CTEST_CONFIGURATION_TYPE}_${CNAME})
endif()

# Empty previous directory before building (otherwise builds are incremental)
if(EMPTY_BINARY_DIRECTORY)
  ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})
endif()

# The Generator for CMake
# ("Unix Makefiles", "Ninja", "Xcode", "NMake Makefiles", "NMake Makefiles JOM",
#  "MinGW Makefiles", "Visual Studio 12 2013", "Visual Studio 14 2015",
#  "Visual Studio 14 2015 Win64", and so on)
if (NOT DEFINED CTEST_CMAKE_GENERATOR)
  file(WRITE ${CTEST_BINARY_DIRECTORY}/cgtest/CMakeLists.txt "message(\"\${CMAKE_GENERATOR}\")\n")
  execute_process(COMMAND ${CMAKE_COMMAND} -Wno-dev .
    ERROR_VARIABLE CTEST_CMAKE_GENERATOR
    OUTPUT_QUIET
    WORKING_DIRECTORY ${CTEST_BINARY_DIRECTORY}/cgtest/
    ERROR_STRIP_TRAILING_WHITESPACE)
endif()

# Cleanup previous tests, settings and test data
file(REMOVE_RECURSE ${CTEST_BINARY_DIRECTORY}/bin)
file(REMOVE ${CTEST_BINARY_DIRECTORY}/CMakeCache.txt)
file(REMOVE_RECURSE ${CTEST_BINARY_DIRECTORY}/Testing)

if (DO_COVERAGE)
  if(NOT DEFINED CTEST_COVERAGE_COMMAND)
    find_program(CTEST_COVERAGE_COMMAND NAMES gcov)
  endif()
    set(CTEST_CUSTOM_COVERAGE_EXCLUDE "${CTEST_SOURCE_DIRECTORY}/external/")
  set(ENV{CXXFLAGS} "$ENV{CXXFLAGS} -g -O0 --coverage -fprofile-arcs -ftest-coverage")
  set(ENV{CFLAGS} "$ENV{CFLAGS} -g -O0 --coverage -fprofile-arcs -ftest-coverage")
endif()

if ("x${CTEST_MEMORYCHECK_TYPE}" STREQUAL "xValgrind")
  find_program(CTEST_MEMORYCHECK_COMMAND NAMES valgrind)
  set(MEMORYCHECK_SUPPRESSIONS_FILE "${CTEST_SOURCE_DIRECTORY}/cmake/valgrind_supp.txt")
  set(MEMORYCHECK_COMMAND_OPTIONS "--error-exitcode=1 --leak-check=yes -q")
  #--tool=memcheck --show-reachable=yes --num-callers=50 --track-origins=yes --trace-children=yes
endif()

if("x${CTEST_MEMORYCHECK_TYPE}" STREQUAL "xAddressSanitizer")
  #See https://github.com/google/sanitizers
  set(ENV{CXXFLAGS} "$ENV{CXXFLAGS} -fsanitize=address -fno-omit-frame-pointer")
  set(ENV{LDFLAGS}  "$ENV{LDFLAGS} -fsanitize=address")
  set(ENV{ASAN_OPTIONS}  "symbolize=1:detect_leaks=1")
  #find_program(LLVMSYM NAMES llvm-symbolizer) #needed in path
  #set(ENV{ASAN_SYMBOLIZER_PATH}  "${LLVMSYM}")
endif()
if("x${CTEST_MEMORYCHECK_TYPE}" STREQUAL "xLeakSanitizer")#part of AddressSanitizer
  set(ENV{CXXFLAGS} "$ENV{CXXFLAGS} -fsanitize=leak -fno-omit-frame-pointer")
  set(ENV{LDFLAGS}  "$ENV{LDFLAGS} -fsanitize=leak")
  #set(ENV{LSAN_OPTIONS} "suppressions="${CTEST_SOURCE_DIRECTORY}/cmake/asan_supp.txt")
endif()
if("x${CTEST_MEMORYCHECK_TYPE}" STREQUAL "xMemorySanitizer")
  # Note: requires full code (including libc++) compiled with -fsanitize=memory
  set(ENV{CXXFLAGS} "$ENV{CXXFLAGS} -fsanitize=memory -fno-omit-frame-pointer -fsanitize-memory-track-origins")
  set(ENV{LDFLAGS}  "$ENV{LDFLAGS} -fsanitize=memory -fsanitize-memory-track-origins")
  set(ENV{MSAN_OPTIONS}  "symbolize=1")
  #find_program(LLVMSYM NAMES llvm-symbolizer) #needed in path
  #set(ENV{MSAN_SYMBOLIZER_PATH}  "${LLVMSYM}")
endif()
if("x${CTEST_MEMORYCHECK_TYPE}" STREQUAL "xThreadSanitizer")
  set(ENV{CXXFLAGS} "$ENV{CXXFLAGS} -fsanitize=thread")
  set(ENV{LDFLAGS}  "$ENV{LDFLAGS} -fsanitize=thread")
  #set(ENV{TSAN_OPTIONS} "report_bugs=1")
endif()
if("x${CTEST_MEMORYCHECK_TYPE}" STREQUAL "xUndefinedBehaviorSanitizer")
  set(ENV{CXXFLAGS} "$ENV{CXXFLAGS} -fsanitize=undefined -fno-omit-frame-pointer")
  set(ENV{LDFLAGS}  "$ENV{LDFLAGS} -fsanitize=undefined")
  set(ENV{UBSAN_OPTIONS} "print_stacktrace=1")
endif()

# Update type (git, svn, wget or url)
if (NOT DEFINED UPDATE_TYPE)
  set(UPDATE_TYPE git)
endif()

if (NOT DEFINED GISMO_BRANCH) #for initial checkout
  set(GISMO_BRANCH stable)
endif()

# Update modules with fetch HEAD commits for all initialized
# submodules
if (NOT DEFINED UPDATE_MODULES)
  set(UPDATE_MODULES OFF)
endif()

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

find_program(CTEST_UPDATE_COMMAND NAMES ${UPDATE_TYPE} ${UPDATE_TYPE}.exe)

# Initial checkout
if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
  if("x${UPDATE_TYPE}" STREQUAL "xgit")
    if("x${UPDATE_PROT}" STREQUAL "xhttps")
      set(gismo_url https://github.com/gismo/gismo.git)
    else() #ssh
      set(gismo_url git@github.com:gismo/gismo.git)
    endif()
    execute_process(COMMAND ${CTEST_UPDATE_COMMAND} clone --depth 1 --branch ${GISMO_BRANCH} ${gismo_url} ${CTEST_SOURCE_DIRECTORY})
    unset(CTEST_CHECKOUT_COMMAND)

  elseif("x${UPDATE_TYPE}" STREQUAL "xsvn")
    if("x${GISMO_BRANCH}" STREQUAL "xstable") # stable
      set(CTEST_CHECKOUT_COMMAND "${CTEST_UPDATE_COMMAND} checkout https://github.com/gismo/gismo.git/trunk ${CTEST_SOURCE_DIRECTORY}")
    else() # branch
      set(CTEST_CHECKOUT_COMMAND "${CTEST_UPDATE_COMMAND} checkout https://github.com/gismo/gismo.git/branches/${GISMO_BRANCH} ${CTEST_SOURCE_DIRECTORY}")
    endif()
  elseif("x${UPDATE_TYPE}" STREQUAL "xwget")
    execute_process(COMMAND /bin/bash "-c" "wget --no-check-certificate -qO - https://github.com/gismo/gismo/archive/${GISMO_BRANCH}.tar.gz | tar -zxf -")
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink gismo-${GISMO_BRANCH} ${CTEST_SOURCE_DIRECTORY})
    set(CTEST_CHECKOUT_COMMAND "${CMAKE_COMMAND} --version")
  elseif("x${UPDATE_TYPE}" STREQUAL "xurl")
    file(DOWNLOAD https://github.com/gismo/gismo/archive/${GISMO_BRANCH}.tar.gz ${CTEST_SCRIPT_DIRECTORY}/${GISMO_BRANCH}.tar.gz)
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E tar xzf ${GISMO_BRANCH}.tar.gz
      COMMAND ${CMAKE_COMMAND} -E create_symlink gismo-${GISMO_BRANCH} ${CTEST_SOURCE_DIRECTORY}
      WORKING_DIRECTORY ${CTEST_SCRIPT_DIRECTORY} )
    set(CTEST_CHECKOUT_COMMAND "${CMAKE_COMMAND} --version")
  endif()
endif()

if("x${UPDATE_TYPE}" STREQUAL "xgit")

  if (NOT "x${GISMO_SUBMODULES}" STREQUAL "x")
    foreach (submod ${GISMO_SUBMODULES})
      #string(TOUPPER ${submod} csubmod)
      #set(SUBM_ARGS ${SUBM_ARGS} -D${csubmod}=ON)
      if ("x${submod}" STREQUAL "xunsupported")
	set(SUBM_ARGS ${SUBM_ARGS} -DGISMO_UNSUPPORTED=ON)
      endif()
      if ("x${submod}" STREQUAL "xmotor")
	set(SUBM_ARGS ${SUBM_ARGS} -DGISMO_MOTOR=ON)
      endif()
      if ("x${submod}" STREQUAL "xgsElasticity")
	set(SUBM_ARGS ${SUBM_ARGS} -DGISMO_ELASTICITY=ON)#GSELASTICITY=ON
      endif()
      if ("x${submod}" STREQUAL "xgsExastencils")
	set(SUBM_ARGS ${SUBM_ARGS} -DGISMO_EXASTENCILS=ON)
      endif()
    endforeach()
  endif()

  foreach (submodule ${GISMO_SUBMODULES})
    if( NOT EXISTS "${CTEST_SOURCE_DIRECTORY}/extensions/${submodule}/.git" )
      execute_process(COMMAND ${CTEST_UPDATE_COMMAND} submodule update --init extensions/${submodule}
	WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY})
    endif()
    if( NOT EXISTS "${CTEST_SOURCE_DIRECTORY}/extensions/${submodule}/.git" )
      message(SEND_ERROR "Problem fetching ${submodule}")
    endif()

    if(${UPDATE_MODULES})
      execute_process(COMMAND ${CTEST_UPDATE_COMMAND} checkout master
	WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY}/extensions/${submodule})
    endif()
  endforeach()
  if(${UPDATE_MODULES})
    set(CTEST_GIT_UPDATE_CUSTOM ${CTEST_UPDATE_COMMAND} pull)
    unset(CTEST_GIT_UPDATE_OPTIONS)
  endif()
endif()

if("${CTEST_CMAKE_GENERATOR}" MATCHES "Make" OR "${CTEST_CMAKE_GENERATOR}" MATCHES "Ninja")
  set(CTEST_USE_LAUNCHERS 1)
else()
  set(CTEST_USE_LAUNCHERS 0)
endif()

set(ENV{CTEST_OUTPUT_ON_FAILURE} 1)
set( $ENV{LC_MESSAGES} "en_EN")
set(ENV{LC_ALL} C)# avoid non-ascii characters

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
#  find_program(UNAME NAMES uname)
#  execute_process(COMMAND "${UNAME}" "-s" OUTPUT_VARIABLE osname OUTPUT_STRIP_TRAILING_WHITESPACE)
#  execute_process(COMMAND "${UNAME}" "-m" OUTPUT_VARIABLE "cpu" OUTPUT_STRIP_TRAILING_WHITESPACE)
#  set(CTEST_BUILD_NAME "${osname}-${cpu} ${CTEST_CMAKE_GENERATOR}-${CTEST_CONFIGURATION_TYPE}-${CNAME}")
  if(${UPDATE_MODULES})
    set(smHead "(head)")
  endif()
  get_filename_component(cxxnamewe "${CXXNAME}" NAME_WE)
  set(CTEST_BUILD_NAME "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR} ${CTEST_CMAKE_GENERATOR}-${CTEST_CONFIGURATION_TYPE}-${cxxnamewe}${smHead}")
  endif()
  message("NAME: ${CTEST_BUILD_NAME}")
  
if(NOT CTEST_BUILD_JOBS)
  include(ProcessorCount)
  ProcessorCount(NPROC)
  #message("Number of processors: ${NPROC}")
  if(${NPROC} EQUAL 0)
    set(NPROC 1)
  endif()
  if(${NPROC} GREATER 40)
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

set(ENV{CTEST_USE_LAUNCHERS_DEFAULT} 1)

macro(get_git_status res)
  if(EXISTS "${CTEST_SOURCE_DIRECTORY}/.git" )
#    execute_process(COMMAND ${CTEST_UPDATE_COMMAND} rev-parse --verify HEAD
#      WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY}
#      OUTPUT_STRIP_TRAILING_WHITESPACE
#      OUTPUT_VARIABLE gitHash)
    execute_process(COMMAND ${CTEST_UPDATE_COMMAND} log -1
      WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY}
      OUTPUT_STRIP_TRAILING_WHITESPACE
      OUTPUT_VARIABLE commitMessage)
    execute_process(COMMAND ${CTEST_UPDATE_COMMAND} submodule status
      WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY}
      OUTPUT_STRIP_TRAILING_WHITESPACE
      OUTPUT_VARIABLE submoduleHashes)
    execute_process(COMMAND ${CTEST_UPDATE_COMMAND} submodule summary
      WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY}
      OUTPUT_STRIP_TRAILING_WHITESPACE
      OUTPUT_VARIABLE submoduleSummary)
    set(${res} "${commitMessage}\n\nSubmodule status:\n${submoduleHashes}\n\nSubmodule summary:\n${submoduleSummary}\n")
  endif()
endmacro(get_git_status)

macro(update_gismo ug_ucount)
  # Reset CTestConfig variables
  if(DEFINED PROJECT_NAME)
    set(CTEST_PROJECT_NAME ${PROJECT_NAME})
    if(NOT DEFINED DROP_LOCATION)
      set(DROP_LOCATION "/submit.php?project=${PROJECT_NAME}")
    endif()
  endif()
  if(DEFINED DROP_LOCATION)
    set(CTEST_DROP_LOCATION ${DROP_LOCATION})
  endif()
  if(DEFINED DROP_SITE)
    set(CTEST_DROP_SITE ${DROP_SITE})
  endif()
  if(DEFINED DROP_METHOD)
    set(CTEST_DROP_METHOD ${DROP_METHOD})
  endif()
  if ("${CMAKE_VERSION}" VERSION_GREATER "3.9.99")
    set(CTEST_LABELS_FOR_SUBPROJECTS ${LABELS_FOR_SUBPROJECTS}) #labels/subprojects
  endif()

  set(ug_updlog "0")
  ctest_update(SOURCE ${CTEST_SOURCE_DIRECTORY} RETURN_VALUE ${ug_ucount})
  ctest_submit(PARTS Update RETRY_COUNT 3 RETRY_DELAY 3)
  set(ug_updlog " ${${ug_ucount}} gismo\n")
  if(${UPDATE_MODULES})
    foreach (submodule ${GISMO_SUBMODULES})
      execute_process(COMMAND ${CTEST_UPDATE_COMMAND} checkout master
	WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY}/extensions/${submodule})
      ctest_update(SOURCE ${CTEST_SOURCE_DIRECTORY}/extensions/${submodule} RETURN_VALUE ug_upd_sm)
      set(ug_updlog "${ug_updlog} ${ug_upd_sm} extensions/${submodule}\n")
      if (${ug_upd_sm} GREATER 0)
	math(EXPR ${ug_ucount} "${${ug_ucount}} + ${ug_upd_sm}")
      endif()
      if(${ug_upd_sm} LESS 0)
        message(SEND_ERROR "Git update submodule error")
      endif()
    endforeach()
  endif()
  get_git_status(gitstatus)
  file(WRITE ${CTEST_BINARY_DIRECTORY}/gitstatus.txt "Revision:\n${gitstatus}\nUpdates:\n${ug_updlog}")
endmacro(update_gismo)

macro(run_ctests)

  # Reset CTestConfig variables
  if(DEFINED PROJECT_NAME)
    set(CTEST_PROJECT_NAME ${PROJECT_NAME})
    if(NOT DEFINED DROP_LOCATION)
      set(DROP_LOCATION "/submit.php?project=${PROJECT_NAME}")
    endif()
  endif()
  if(DEFINED DROP_LOCATION)
    set(CTEST_DROP_LOCATION ${DROP_LOCATION})
  endif()
  if(DEFINED DROP_SITE)
    set(CTEST_DROP_SITE ${DROP_SITE})
  endif()
  if(DEFINED DROP_METHOD)
    set(CTEST_DROP_METHOD ${DROP_METHOD})
  endif()
  if ("${CMAKE_VERSION}" VERSION_GREATER "3.9.99")
    set(CTEST_LABELS_FOR_SUBPROJECTS ${LABELS_FOR_SUBPROJECTS}) #labels/subprojects
  endif()
  
  ctest_configure(OPTIONS "${CMAKE_ARGS};${SUBM_ARGS};-DCTEST_USE_LAUNCHERS=${CTEST_USE_LAUNCHERS};-DBUILD_TESTING=ON;-DDART_TESTING_TIMEOUT=${CTEST_TEST_TIMEOUT}")
  #ctest_submit(PARTS Configure Update  RETRY_COUNT 3 RETRY_DELAY 3)
  ctest_submit(PARTS Configure Notes RETRY_COUNT 3 RETRY_DELAY 3)

  #"${CMAKE_VERSION}" VERSION_LESS "3.10"
  if(NOT "x${LABELS_FOR_SUBPROJECTS}" STREQUAL "x")

    foreach(subproject ${LABELS_FOR_SUBPROJECTS})
      message("Subproject ${subproject}")
      if ("${CMAKE_VERSION}" VERSION_LESS "3.10")
        set_property(GLOBAL PROPERTY SubProject ${subproject})
        set_property(GLOBAL PROPERTY Label ${subproject})
      endif()
      ctest_build(TARGET ${subproject} APPEND)
      ctest_submit(PARTS Build  RETRY_COUNT 3 RETRY_DELAY 3)
      ctest_test(INCLUDE_LABEL "${subproject}" PARALLEL_LEVEL ${CTEST_TEST_JOBS})
      ctest_submit(PARTS Test  RETRY_COUNT 3 RETRY_DELAY 3)

      if(DO_COVERAGE)
        ctest_coverage(BUILD "${CTEST_BINARY_DIRECTORY}" LABELS "${subproject}" APPEND)
        ctest_submit(PARTS Coverage  RETRY_COUNT 3 RETRY_DELAY 3)
      endif()

      if(NOT "x${CTEST_MEMORYCHECK_TYPE}" STREQUAL "xNone")
        ctest_memcheck(INCLUDE_LABEL "${subproject}" APPEND)
        ctest_submit(PARTS MemCheck  RETRY_COUNT 3 RETRY_DELAY 3)
      endif()

    endforeach()

  else() # No subprojects


    if("x${CTEST_CMAKE_GENERATOR}" STREQUAL "xNinja")
      ctest_build(TARGET UnitTestPP APPEND) # for older versions of ninja
    endif()
    ctest_submit(PARTS Build  RETRY_COUNT 3 RETRY_DELAY 3)
    ctest_build(APPEND)
    ctest_submit(PARTS Build  RETRY_COUNT 3 RETRY_DELAY 3)
    ctest_build(TARGET unittests APPEND)
    ctest_submit(PARTS Build  RETRY_COUNT 3 RETRY_DELAY 3)
    ctest_test(PARALLEL_LEVEL ${CTEST_TEST_JOBS})
    ctest_submit(PARTS Test  RETRY_COUNT 3 RETRY_DELAY 3)

    if(DO_COVERAGE)
      #message("Running coverage..")
      ctest_coverage(BUILD "${CTEST_BINARY_DIRECTORY}" APPEND)
      ctest_submit(PARTS Coverage  RETRY_COUNT 3 RETRY_DELAY 3)
    endif()

    if(NOT "x${CTEST_MEMORYCHECK_TYPE}" STREQUAL "xNone")
      #message("Running memcheck..")
      ctest_memcheck()
      ctest_submit(PARTS MemCheck  RETRY_COUNT 3 RETRY_DELAY 3)
    endif()

  endif()
endmacro(run_ctests)

file(MAKE_DIRECTORY "${CTEST_BINARY_DIRECTORY}")

set(CTEST_NOTES_FILES ${CTEST_BINARY_DIRECTORY}/gitstatus.txt)

if(NOT "${CTEST_TEST_MODEL}" STREQUAL "Continuous")

  ctest_start(${CTEST_TEST_MODEL})
  if(NOT "${CTEST_UPDATE_COMMAND}" STREQUAL "CTEST_UPDATE_COMMAND-NOTFOUND")
    update_gismo(updcount)
  endif()
  run_ctests()

else() #continuous model

  while(${CTEST_ELAPSED_TIME} LESS ${test_runtime})
    set(START_TIME ${CTEST_ELAPSED_TIME})
    ctest_start(${CTEST_TEST_MODEL})
    update_gismo(updcount)
    if( ${updcount} GREATER 0 )
      run_ctests()
    endif()
    ctest_sleep(${START_TIME} 300 ${CTEST_ELAPSED_TIME})
  endwhile()

endif(NOT "${CTEST_TEST_MODEL}" STREQUAL "Continuous")
