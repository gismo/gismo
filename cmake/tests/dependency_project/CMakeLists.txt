######################################################################
## dependency_project/CMakeLists.txt
##
## The project under test for run_dependency_tests.cmake: sandboxes
## gismo_SOURCE_DIR/gismo_BINARY_DIR, includes the real gsFetch.cmake, makes
## exactly one gismo_add_dependency() call selected by -DCASE=<id>, and
## records the observables the driver asserts on into dep_report.txt. A
## NEVER-mode case never reaches the report: every failure path inside
## gismo_add_dependency() is message(FATAL_ERROR ...), which aborts this
## configure outright - those cases are asserted by the driver directly on
## the sub-configure's non-zero RESULT_VARIABLE and captured stdout+stderr.
######################################################################

if(CMAKE_VERSION VERSION_LESS "3.19")
  cmake_minimum_required(VERSION 2.8.12)
else()
  cmake_minimum_required(VERSION 3.1...3.10)
endif()
project(gsdeptest CXX)

foreach(_gsdep_required CASE SANDBOX_SRC_DIR SANDBOX_BIN_DIR GISMO_CMAKE_DIR)
  if(NOT DEFINED ${_gsdep_required})
    message(FATAL_ERROR "gsdeptest: -D${_gsdep_required}=... is required")
  endif()
endforeach()

## Sandbox both globals gismo_fetch_directory() reads/writes. A function()
## reads its caller's scope, so setting these plain (non-cache) variables
## here - before gsFetch.cmake is include()d and before any
## gismo_add_dependency() call - is what keeps every fetch/vendor
## destination and every SEND_ERROR path confined to the scratch tree; see
## gsFetch.cmake:26-38.
set(gismo_SOURCE_DIR "${SANDBOX_SRC_DIR}")
set(gismo_BINARY_DIR "${SANDBOX_BIN_DIR}")

include("${GISMO_CMAKE_DIR}/gsFetch.cmake")

## Normalises a semicolon-separated list into the report's '|'-joined form,
## realpath'ing each plain-path element individually (a still-unresolved
## generator expression, "$<...>", is kept verbatim - REALPATH on literal
## genex text is meaningless). Used for GISMO_INCLUDE_DIRS both before and
## after the gismo_add_dependency() call, so the driver's own REALPATH of
## its expected value compares like for like even when SCRATCH_DIR sits
## behind a symlink (e.g. macOS /tmp -> /private/tmp).
function(gsdep_format_dir_list VALUE OUT)
  if(NOT DEFINED VALUE)
    set(${OUT} "<unset>" PARENT_SCOPE)
    return()
  endif()
  set(_norm "")
  foreach(_d ${VALUE})
    if(_d MATCHES "\\$<")
      list(APPEND _norm "${_d}")
    else()
      get_filename_component(_d_real "${_d}" REALPATH)
      list(APPEND _norm "${_d_real}")
    endif()
  endforeach()
  string(REPLACE ";" "|" _joined "${_norm}")
  set(${OUT} "${_joined}" PARENT_SCOPE)
endfunction()

if(DEFINED GISMO_INCLUDE_DIRS)
  gsdep_format_dir_list("${GISMO_INCLUDE_DIRS}" _gsdep_before_dirs)
else()
  set(_gsdep_before_dirs "<unset>")
endif()
if(DEFINED gismo_LINKER)
  string(REPLACE ";" "|" _gsdep_before_linker "${gismo_LINKER}")
else()
  set(_gsdep_before_linker "<unset>")
endif()

set(_gsdep_name "")
set(_gsdep_consumer "")

if(CASE STREQUAL "header_found")
  set(_gsdep_name FakeHdrDep)
  gismo_add_dependency(FakeHdrDep
    FIND_PACKAGE FakeHdrDep
    MODE HEADER_ONLY
    TARGET GsDepTest::FakeHdrDep)
  set(_gsdep_consumer use_header)

elseif(CASE STREQUAL "header_vendored")
  set(_gsdep_name VendoredHdrDep)
  gismo_add_dependency(VendoredHdrDep
    GIT_REPOSITORY https://invalid.invalid/gsdeptest-should-never-be-cloned.git
    MODE HEADER_ONLY
    INCLUDE_SUBDIR include
    TARGET GsDepTest::VendoredHdrDep)
  set(_gsdep_consumer use_header)

elseif(CASE STREQUAL "header_never")
  set(_gsdep_name MissingHdrDep)
  gismo_add_dependency(MissingHdrDep
    MODE HEADER_ONLY
    TARGET GsDepTest::MissingHdrDep)

elseif(CASE STREQUAL "sources_found")
  set(_gsdep_name FakeSrcDep)
  gismo_add_dependency(FakeSrcDep
    FIND_PACKAGE FakeSrcDep
    MODE SOURCES
    SOURCE_GLOBS "src/*.cpp"
    TARGET GsDepTest::FakeSrcDep)
  set(_gsdep_consumer use_sources)

elseif(CASE STREQUAL "sources_vendored")
  set(_gsdep_name VendoredSrcDep)
  gismo_add_dependency(VendoredSrcDep
    GIT_REPOSITORY https://invalid.invalid/gsdeptest-should-never-be-cloned.git
    MODE SOURCES
    INCLUDE_SUBDIR include
    SOURCE_GLOBS "src/*.cpp"
    TARGET GsDepTest::VendoredSrcDep)
  set(_gsdep_consumer use_sources)

elseif(CASE STREQUAL "sources_never")
  set(_gsdep_name MissingSrcDep)
  gismo_add_dependency(MissingSrcDep
    MODE SOURCES
    SOURCE_GLOBS "src/*.cpp"
    TARGET GsDepTest::MissingSrcDep)

elseif(CASE STREQUAL "header_network")
  # No INCLUDE_SUBDIR: gsUnitTest keeps its headers at the repository root.
  set(_gsdep_name NetHdrDep)
  gismo_add_dependency(NetHdrDep
    GIT_REPOSITORY https://github.com/gismo/gsUnitTest.git
    MODE HEADER_ONLY
    TARGET GsDepTest::NetHdrDep)

elseif(CASE STREQUAL "found_multi_twice")
  set(_gsdep_name MultiTwiceDep)
  gismo_add_dependency(MultiTwiceDep
    FIND_PACKAGE MultiTwiceDep
    MODE HEADER_ONLY
    TARGET GsDepTest::MultiTwiceDep)

elseif(CASE STREQUAL "found_no_include")
  set(_gsdep_name NoIncludeDep)
  gismo_add_dependency(NoIncludeDep
    FIND_PACKAGE NoIncludeDep
    MODE HEADER_ONLY
    TARGET GsDepTest::NoIncludeDep)

elseif(CASE STREQUAL "genex_build_multi")
  set(_gsdep_name GenexBuildMultiDep)
  gismo_add_dependency(GenexBuildMultiDep
    FIND_PACKAGE GenexBuildMultiDep
    MODE HEADER_ONLY
    TARGET GsDepTest::GenexBuildMultiDep)

elseif(CASE STREQUAL "genex_concat")
  set(_gsdep_name GenexConcatDep)
  gismo_add_dependency(GenexConcatDep
    FIND_PACKAGE GenexConcatDep
    MODE HEADER_ONLY
    TARGET GsDepTest::GenexConcatDep)

elseif(CASE STREQUAL "genex_install_only")
  set(_gsdep_name GenexInstallOnlyDep)
  gismo_add_dependency(GenexInstallOnlyDep
    FIND_PACKAGE GenexInstallOnlyDep
    MODE HEADER_ONLY
    TARGET GsDepTest::GenexInstallOnlyDep)

elseif(CASE STREQUAL "genex_unknown")
  set(_gsdep_name GenexUnknownDep)
  gismo_add_dependency(GenexUnknownDep
    FIND_PACKAGE GenexUnknownDep
    MODE HEADER_ONLY
    TARGET GsDepTest::GenexUnknownDep)

else()
  message(FATAL_ERROR "gsdeptest: unknown CASE '${CASE}'")
endif()

## Reaching this point means gismo_add_dependency() resolved the dependency
## (it never returns on a failure path - see the comment above).
set(_gsdep_target "GsDepTest::${_gsdep_name}")
set(_gsdep_report "${CMAKE_BINARY_DIR}/dep_report.txt")

if(DEFINED GISMO_INCLUDE_DIRS)
  gsdep_format_dir_list("${GISMO_INCLUDE_DIRS}" _gsdep_after_dirs)
else()
  set(_gsdep_after_dirs "<unset>")
endif()
if(DEFINED gismo_LINKER)
  string(REPLACE ";" "|" _gsdep_after_linker "${gismo_LINKER}")
else()
  set(_gsdep_after_linker "<unset>")
endif()

if(TARGET ${_gsdep_target})
  set(_gsdep_target_exists "TRUE")

  get_target_property(_gsdep_target_type ${_gsdep_target} TYPE)
  if(NOT _gsdep_target_type)
    set(_gsdep_target_type "<unset>")
  endif()

  get_target_property(_gsdep_target_include ${_gsdep_target} INTERFACE_INCLUDE_DIRECTORIES)
  if(_gsdep_target_include STREQUAL "_gsdep_target_include-NOTFOUND" OR NOT _gsdep_target_include)
    set(_gsdep_target_include "<unset>")
  else()
    gsdep_format_dir_list("${_gsdep_target_include}" _gsdep_target_include)
  endif()
else()
  set(_gsdep_target_exists "FALSE")
  set(_gsdep_target_type "<unset>")
  set(_gsdep_target_include "<unset>")
endif()

if(DEFINED ${_gsdep_name}_FOUND)
  set(_gsdep_found_val "${${_gsdep_name}_FOUND}")
else()
  set(_gsdep_found_val "<unset>")
endif()
if(DEFINED ${_gsdep_name}_VENDORED)
  set(_gsdep_vendored_val "${${_gsdep_name}_VENDORED}")
else()
  set(_gsdep_vendored_val "<unset>")
endif()

file(WRITE "${_gsdep_report}"
"CASE=${CASE}
DEP_NAME=${_gsdep_name}
TARGET_EXISTS=${_gsdep_target_exists}
TARGET_TYPE=${_gsdep_target_type}
TARGET_INCLUDE=${_gsdep_target_include}
FOUND=${_gsdep_found_val}
VENDORED=${_gsdep_vendored_val}
INCLUDE_DIRS_BEFORE=${_gsdep_before_dirs}
INCLUDE_DIRS_AFTER=${_gsdep_after_dirs}
LINKER_BEFORE=${_gsdep_before_linker}
LINKER_AFTER=${_gsdep_after_linker}
")

if(NOT _gsdep_consumer STREQUAL "")
  if(NOT DEFINED GSDEP_CONSUMER_DIR)
    message(FATAL_ERROR "gsdeptest: -DGSDEP_CONSUMER_DIR=... is required for CASE=${CASE}")
  endif()
  # Pin the runtime output directory to the same path for every
  # configuration, single- or multi-config generator alike, so the driver
  # never has to guess a per-config subdirectory when it looks for the
  # built executable.
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
  foreach(_gsdep_cfg Debug Release RelWithDebInfo MinSizeRel)
    string(TOUPPER "${_gsdep_cfg}" _gsdep_cfg_upper)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${_gsdep_cfg_upper} "${CMAKE_BINARY_DIR}/bin")
  endforeach()
  add_executable(${_gsdep_consumer} "${GSDEP_CONSUMER_DIR}/${_gsdep_consumer}.cpp")
  target_link_libraries(${_gsdep_consumer} PRIVATE ${_gsdep_target})
endif()
