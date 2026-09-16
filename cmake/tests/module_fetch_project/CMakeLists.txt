######################################################################
## module_fetch_project/CMakeLists.txt
##
## The project under test for run_module_fetch_tests.cmake: sandboxes
## gismo_SOURCE_DIR/gismo_BINARY_DIR, includes the real gismoUse.cmake (for
## get_repo_info(), used by gismo_fetch_module_source()) and gsFetch.cmake,
## sets GISMO_OPTIONAL to a single fixture module ("gsFetchFixture") fetched
## from a local file:// git repository (no network), runs exactly one
## -DCASE=<id> scenario, and records the observables the driver asserts on
## into module_fetch_report.txt.
##
## Nothing is compiled here (project(... NONE)) - gismo_add_extension() only
## add_subdirectory()s the fixture's own CMakeLists.txt, which defines a
## single add_custom_target().
######################################################################

cmake_minimum_required(VERSION 3.1...3.10)
project(gsmodfetchtest NONE)

foreach(_gsmf_required CASE SANDBOX_SRC_DIR SANDBOX_BIN_DIR GISMO_CMAKE_DIR gsFetchFixture_url)
  if(NOT DEFINED ${_gsmf_required})
    message(FATAL_ERROR "gsmodfetchtest: -D${_gsmf_required}=... is required")
  endif()
endforeach()

## gismo_fetch_module_source()'s clone branch requires GIT_EXECUTABLE; the
## present-module update branch (not exercised here, since the fixture is
## always absent before the first fetch) also reads it directly without a
## find_package() call of its own.
find_package(Git REQUIRED)

## Sandbox both globals gismo_fetch_module_source()/gismo_add_extension()
## read/write - see cmake/tests/dependency_project/CMakeLists.txt for
## the same rationale.
set(gismo_SOURCE_DIR "${SANDBOX_SRC_DIR}")
set(gismo_BINARY_DIR "${SANDBOX_BIN_DIR}")

## The clone in gismo_fetch_module_source() runs with
## WORKING_DIRECTORY ${gismo_SOURCE_DIR}/optional, which must already exist.
file(MAKE_DIRECTORY "${gismo_SOURCE_DIR}/optional")

include("${GISMO_CMAKE_DIR}/gismoUse.cmake")
include("${GISMO_CMAKE_DIR}/gsFetch.cmake")

set(GISMO_OPTIONAL gsFetchFixture)
## Every captured include-dir list below is non-empty from the start, so
## "the fixture include dir is absent" can never be satisfied vacuously by a
## list that was never populated.
set(GISMO_INCLUDE_DIRS "${SANDBOX_SRC_DIR}/sentinel_include")

## Same '|'-joined, REALPATH'd report formatting as
## cmake/tests/dependency_project/CMakeLists.txt.
function(gsmf_format_dir_list VALUE OUT)
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

set(_gsmf_exists_before FALSE)
if(EXISTS "${gismo_SOURCE_DIR}/optional/gsFetchFixture/CMakeLists.txt")
  set(_gsmf_exists_before TRUE)
endif()

## Defaults for the two case-specific keys; overwritten in the branches that
## produce them so every case's report carries every key.
set(_gsmf_include_dirs_after_hook "<unset>")
set(_gsmf_marker_survives "<unset>")

## Case dispatch: deliberately top-level if/elseif, never wrapped in a
## function() - gismo_prepare_optional_modules()/gismo_include_module_dependencies()
## are macros and dependencies.cmake sets plain variables (GISMO_INCLUDE_DIRS),
## which a function() scope would swallow.
if(CASE STREQUAL "absent_module")
  gismo_prepare_optional_modules()
  gsmf_format_dir_list("${GISMO_INCLUDE_DIRS}" _gsmf_include_dirs_final)

elseif(CASE STREQUAL "legacy_order")
  ## The legacy top-level loop: include a listed module's dependencies.cmake
  ## if it already exists, then fetch the module. Replayed here with the one
  ## substitution PROJECT_SOURCE_DIR -> gismo_SOURCE_DIR required to keep the
  ## falsifier non-vacuous in this sandbox (PROJECT_SOURCE_DIR here is
  ## cmake/tests/module_fetch_project, where optional/gsFetchFixture never
  ## exists).
  list(REMOVE_DUPLICATES GISMO_OPTIONAL)
  list(SORT GISMO_OPTIONAL)
  foreach(gssm ${GISMO_OPTIONAL})
    string(STRIP ${gssm} SUBMODULE)
    if(EXISTS "${gismo_SOURCE_DIR}/optional/${SUBMODULE}/dependencies.cmake")
      message(STATUS "Processing dependencies for optional module: ${SUBMODULE}")
      include("${gismo_SOURCE_DIR}/optional/${SUBMODULE}/dependencies.cmake")
    endif()
  endforeach()
  gsmf_format_dir_list("${GISMO_INCLUDE_DIRS}" _gsmf_include_dirs_after_hook)
  gismo_fetch_module(gsFetchFixture)
  gsmf_format_dir_list("${GISMO_INCLUDE_DIRS}" _gsmf_include_dirs_final)

elseif(CASE STREQUAL "idempotent")
  ## GISMO_SUBMODULES_HEAD=ON makes the present-module branch of a second,
  ## unguarded gismo_fetch_module_source() call run
  ## fetch/reset --hard/clean -dfx/rebase, which would delete the untracked
  ## marker below - so MARKER_SURVIVES actually discriminates a missing
  ## fetch guard instead of surviving trivially.
  gismo_prepare_optional_modules()
  file(WRITE "${gismo_SOURCE_DIR}/optional/gsFetchFixture/untracked_marker.txt" "x")
  gismo_prepare_optional_modules()
  gismo_fetch_module(gsFetchFixture)
  gsmf_format_dir_list("${GISMO_INCLUDE_DIRS}" _gsmf_include_dirs_final)
  set(_gsmf_marker_survives FALSE)
  if(EXISTS "${gismo_SOURCE_DIR}/optional/gsFetchFixture/untracked_marker.txt")
    set(_gsmf_marker_survives TRUE)
  endif()

else()
  message(FATAL_ERROR "gsmodfetchtest: unknown CASE '${CASE}'")
endif()

set(_gsmf_cloned FALSE)
if(EXISTS "${gismo_SOURCE_DIR}/optional/gsFetchFixture/CMakeLists.txt")
  set(_gsmf_cloned TRUE)
endif()

get_property(_gsmf_deps_included GLOBAL PROPERTY GSFETCHFIXTURE_DEPS_INCLUDED)
if(NOT _gsmf_deps_included)
  set(_gsmf_deps_included 0)
endif()

file(WRITE "${CMAKE_BINARY_DIR}/module_fetch_report.txt"
"CASE=${CASE}
EXISTS_BEFORE=${_gsmf_exists_before}
CLONED=${_gsmf_cloned}
INCLUDE_DIRS_AFTER_HOOK=${_gsmf_include_dirs_after_hook}
INCLUDE_DIRS_FINAL=${_gsmf_include_dirs_final}
DEPS_INCLUDED=${_gsmf_deps_included}
MARKER_SURVIVES=${_gsmf_marker_survives}
")
