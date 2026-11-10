######################################################################
## run_module_fetch_tests.cmake
##
## Standalone driver (invoked as `cmake -P run_module_fetch_tests.cmake`,
## wired into ctest as the single test `cmake_module_fetch_test`) for
## gismo_prepare_optional_modules()/gismo_fetch_module_source()/
## gismo_include_module_dependencies() (cmake/gsFetch.cmake). Builds a
## fixture git repository once under ${SCRATCH_DIR} and runs three
## sub-configures of cmake/tests/module_fetch_project against it:
##
##  - absent_module: gismo_prepare_optional_modules() fetches a module
##    missing from optional/ BEFORE its dependencies.cmake is read, so the
##    module's include dir reaches GISMO_INCLUDE_DIRS.
##  - legacy_order: the falsifier - replaying the legacy top-level order
##    (include the module's dependencies.cmake if it exists, then
##    gismo_fetch_module()) on the same fixture does NOT pick up the include
##    dir, because the hook runs before the module exists.
##  - idempotent: two consecutive gismo_prepare_optional_modules() calls (and
##    a further gismo_fetch_module()) neither re-clone nor re-include.
##
## No network: the fixture module is cloned from a local file:// git
## repository this driver creates at test time. If no git executable is
## found, the whole ctest is reported SKIPPED (never PASSED) via the
## GSMODFETCH_SKIPPED marker below and SKIP_REGULAR_EXPRESSION at the
## registration site.
######################################################################

cmake_minimum_required(VERSION 3.1...3.10)

foreach(_gsmf_required GISMO_SOURCE_DIR SCRATCH_DIR TEST_GENERATOR TEST_MAKE_PROGRAM)
  if(NOT DEFINED ${_gsmf_required})
    message(FATAL_ERROR "run_module_fetch_tests.cmake: -D${_gsmf_required}=... is required")
  endif()
endforeach()

## SKIP_REGULAR_EXPRESSION (used by the registration to report this ctest as
## SKIPPED rather than PASSED/FAILED) needs CMake >= 3.16. Below that, a
## missing git executable is a loud configure-time failure instead of a
## silent pass.
find_program(_gsmf_git NAMES git git.exe git.cmd)
if(NOT _gsmf_git)
  if(CMAKE_VERSION VERSION_LESS "3.16")
    message(FATAL_ERROR "run_module_fetch_tests.cmake: git executable not found and CMake < 3.16 cannot SKIP this test")
  endif()
  message(STATUS "GSMODFETCH_SKIPPED: git executable not found - cmake_module_fetch_test cannot run")
  return()
endif()

set(GSMF_PROJECT_DIR "${GISMO_SOURCE_DIR}/cmake/tests/module_fetch_project")

file(REMOVE_RECURSE "${SCRATCH_DIR}")
file(MAKE_DIRECTORY "${SCRATCH_DIR}")

set_property(GLOBAL PROPERTY gsmf_fail_count 0)
set_property(GLOBAL PROPERTY gsmf_pass_count 0)
set_property(GLOBAL PROPERTY gsmf_skip_count 0)
set_property(GLOBAL PROPERTY gsmf_case_start_fail 0)
## gsmf_fail_count counts failing *assertions*; gsmf_case_fail_count counts
## failing *cases* and is what the final RESULT: line and the FATAL_ERROR
## gate use - same split as run_dependency_tests.cmake.
set_property(GLOBAL PROPERTY gsmf_case_fail_count 0)

######################################################################
## expectation helpers - each prints exactly
##   FAIL <case>/<label>: expected '<expected>' got '<actual>'
## on mismatch and increments the shared fail counter (a GLOBAL PROPERTY,
## since these run inside function() scopes).
######################################################################

function(_gsmf_inc PROP)
  get_property(_v GLOBAL PROPERTY ${PROP})
  math(EXPR _v "${_v}+1")
  set_property(GLOBAL PROPERTY ${PROP} ${_v})
endfunction()

function(expect_str CASE_ID LABEL EXPECTED ACTUAL)
  if(NOT "${EXPECTED}" STREQUAL "${ACTUAL}")
    message(STATUS "FAIL ${CASE_ID}/${LABEL}: expected '${EXPECTED}' got '${ACTUAL}'")
    _gsmf_inc(gsmf_fail_count)
  endif()
endfunction()

function(expect_bool CASE_ID LABEL EXPECTED_BOOL ACTUAL)
  ## EXPECTED_BOOL is the literal 'TRUE' or 'FALSE'; the report always
  ## writes CMake's own TRUE/FALSE for these keys, but '<unset>' (a key
  ## missing from the report) is always a FAIL - it is never a silent
  ## substitute for FALSE.
  if("${ACTUAL}" STREQUAL "<unset>")
    message(STATUS "FAIL ${CASE_ID}/${LABEL}: expected '${EXPECTED_BOOL}' got '<unset>'")
    _gsmf_inc(gsmf_fail_count)
    return()
  endif()
  if(NOT "${ACTUAL}" STREQUAL "${EXPECTED_BOOL}")
    message(STATUS "FAIL ${CASE_ID}/${LABEL}: expected '${EXPECTED_BOOL}' got '${ACTUAL}'")
    _gsmf_inc(gsmf_fail_count)
  endif()
endfunction()

function(expect_list_contains CASE_ID LABEL EXPECTED_DIR ACTUAL_LIST)
  ## ACTUAL_LIST is a '|'-joined, REALPATH'd report value - see
  ## expect_list_contains in run_dependency_tests.cmake (same shape; the two
  ## drivers are separate scripts).
  if("${ACTUAL_LIST}" STREQUAL "<unset>")
    message(STATUS "FAIL ${CASE_ID}/${LABEL}: expected to contain '${EXPECTED_DIR}' got '<unset>'")
    _gsmf_inc(gsmf_fail_count)
    return()
  endif()
  string(REPLACE "|" ";" _elems "${ACTUAL_LIST}")
  set(_found FALSE)
  foreach(_e ${_elems})
    if("${_e}" STREQUAL "${EXPECTED_DIR}")
      set(_found TRUE)
    endif()
  endforeach()
  if(NOT _found)
    message(STATUS "FAIL ${CASE_ID}/${LABEL}: expected '${EXPECTED_DIR}' got '${ACTUAL_LIST}'")
    _gsmf_inc(gsmf_fail_count)
  endif()
endfunction()

function(expect_list_not_contains CASE_ID LABEL FORBIDDEN_DIR ACTUAL_LIST)
  ## Inverse of expect_list_contains(): '<unset>' is a FAIL here too, since
  ## an unpopulated list can never demonstrate the falsifier's actual claim
  ## (the hook ran before the module existed) - it would merely be a list
  ## that was never captured.
  if("${ACTUAL_LIST}" STREQUAL "<unset>")
    message(STATUS "FAIL ${CASE_ID}/${LABEL}: expected to not contain '${FORBIDDEN_DIR}' got '<unset>'")
    _gsmf_inc(gsmf_fail_count)
    return()
  endif()
  string(REPLACE "|" ";" _elems "${ACTUAL_LIST}")
  foreach(_e ${_elems})
    if("${_e}" STREQUAL "${FORBIDDEN_DIR}")
      message(STATUS "FAIL ${CASE_ID}/${LABEL}: expected to not contain '${FORBIDDEN_DIR}' got '${ACTUAL_LIST}'")
      _gsmf_inc(gsmf_fail_count)
      return()
    endif()
  endforeach()
endfunction()

function(expect_count_matches CASE_ID LABEL PATTERN EXPECTED_COUNT ACTUAL_TEXT)
  string(REGEX MATCHALL "${PATTERN}" _matches "${ACTUAL_TEXT}")
  list(LENGTH _matches _n)
  if(NOT _n EQUAL EXPECTED_COUNT)
    message(STATUS "FAIL ${CASE_ID}/${LABEL}: expected '${EXPECTED_COUNT}' occurrence(s) of '${PATTERN}' got '${_n}'")
    _gsmf_inc(gsmf_fail_count)
  endif()
endfunction()

######################################################################
## case bookkeeping and report parsing - same shape as
## run_dependency_tests.cmake.
######################################################################

function(gsmf_case_begin CASE_ID)
  get_property(_f GLOBAL PROPERTY gsmf_fail_count)
  set_property(GLOBAL PROPERTY gsmf_case_start_fail ${_f})
endfunction()

function(gsmf_case_end CASE_ID)
  get_property(_f0 GLOBAL PROPERTY gsmf_case_start_fail)
  get_property(_f1 GLOBAL PROPERTY gsmf_fail_count)
  if(_f1 GREATER _f0)
    message(STATUS "CASE ${CASE_ID} FAIL")
    _gsmf_inc(gsmf_case_fail_count)
  else()
    message(STATUS "CASE ${CASE_ID} PASS")
    _gsmf_inc(gsmf_pass_count)
  endif()
endfunction()

function(gsmf_read_report FILE OUT_PREFIX)
  file(STRINGS "${FILE}" _lines)
  foreach(_line ${_lines})
    if(_line MATCHES "^([A-Z_]+)=(.*)$")
      set(${OUT_PREFIX}_${CMAKE_MATCH_1} "${CMAKE_MATCH_2}" PARENT_SCOPE)
    endif()
  endforeach()
endfunction()

function(gsmf_check_containment CASE_ID)
  if(EXISTS "${GISMO_SOURCE_DIR}/optional/gsFetchFixture")
    message(STATUS "FAIL ${CASE_ID}/containment: expected 'NOT EXISTS' got 'EXISTS (${GISMO_SOURCE_DIR}/optional/gsFetchFixture)'")
    _gsmf_inc(gsmf_fail_count)
  endif()
endfunction()

## Configures CASE_SRC_DIR into CASE_BIN_DIR, forwarding the generator/make
## program pair plus every extra -D flag given as ARGN. No
## CMAKE_CXX_COMPILER/CMAKE_BUILD_TYPE forwarding (unlike
## run_dependency_tests.cmake's gsdep_configure) - module_fetch_project uses
## project(... NONE), so no compiler is ever detected.
function(gsmf_configure CASE_ID CASE_SRC_DIR CASE_BIN_DIR OUT_RESULT_VAR OUT_OUTPUT_VAR)
  file(MAKE_DIRECTORY "${CASE_BIN_DIR}")
  execute_process(
    COMMAND ${CMAKE_COMMAND}
      -G "${TEST_GENERATOR}"
      -DCMAKE_MAKE_PROGRAM=${TEST_MAKE_PROGRAM}
      ${ARGN}
      "${CASE_SRC_DIR}"
    WORKING_DIRECTORY "${CASE_BIN_DIR}"
    RESULT_VARIABLE _result
    OUTPUT_VARIABLE _output
    ERROR_VARIABLE _output)
  set(${OUT_RESULT_VAR} "${_result}" PARENT_SCOPE)
  set(${OUT_OUTPUT_VAR} "${_output}" PARENT_SCOPE)
endfunction()

######################################################################
## fixture repository - created once under ${SCRATCH_DIR}, shared read-only
## by every case (a clone never writes to its origin). The only git writes
## this driver performs; see the "Git safety" note on each command below.
######################################################################

set(GSMF_FIXTURE_DIR "${SCRATCH_DIR}/fixture_repo/gsFetchFixture")
file(MAKE_DIRECTORY "${GSMF_FIXTURE_DIR}/include")

file(WRITE "${GSMF_FIXTURE_DIR}/CMakeLists.txt" "add_custom_target(gsFetchFixture)\n")

## Increments a GLOBAL property every time this file is include()d (so the
## driver can assert "included exactly once per configure") and appends the
## CLONE's own include dir - CMAKE_CURRENT_LIST_DIR resolves to the
## directory of the file actually being processed, i.e. the clone under
## <sandbox-src>/optional/gsFetchFixture, never this fixture repository.
file(WRITE "${GSMF_FIXTURE_DIR}/dependencies.cmake"
"get_property(_gsfixture_n GLOBAL PROPERTY GSFETCHFIXTURE_DEPS_INCLUDED)
if(NOT _gsfixture_n)
  set(_gsfixture_n 0)
endif()
math(EXPR _gsfixture_n \"\${_gsfixture_n}+1\")
set_property(GLOBAL PROPERTY GSFETCHFIXTURE_DEPS_INCLUDED \${_gsfixture_n})
set(GISMO_INCLUDE_DIRS \${GISMO_INCLUDE_DIRS} \"\${CMAKE_CURRENT_LIST_DIR}/include\")
unset(_gsfixture_n)
")

file(WRITE "${GSMF_FIXTURE_DIR}/include/gsFetchFixture.h" "// gsFetchFixture module-fetch test fixture header.\n")

## Git safety: init is confined to GSMF_FIXTURE_DIR via WORKING_DIRECTORY and
## verified by both a zero result AND an actual .git directory, so a failed
## init can never let add/commit fall through onto the enclosing worktree.
execute_process(COMMAND "${_gsmf_git}" init
  WORKING_DIRECTORY "${GSMF_FIXTURE_DIR}"
  RESULT_VARIABLE _gsmf_git_init_res
  OUTPUT_QUIET ERROR_QUIET)
if(NOT _gsmf_git_init_res EQUAL 0 OR NOT EXISTS "${GSMF_FIXTURE_DIR}/.git")
  message(FATAL_ERROR "run_module_fetch_tests.cmake: git init failed in fixture dir (result ${_gsmf_git_init_res})")
endif()

## add/commit name the repository explicitly via --git-dir/--work-tree so
## they can only ever touch GSMF_FIXTURE_DIR, regardless of the process's
## working directory.
execute_process(COMMAND "${_gsmf_git}"
    "--git-dir=${GSMF_FIXTURE_DIR}/.git" "--work-tree=${GSMF_FIXTURE_DIR}"
    -c user.name=gismo-test -c user.email=test@invalid
    add -A
  RESULT_VARIABLE _gsmf_git_add_res
  OUTPUT_QUIET ERROR_QUIET)
if(NOT _gsmf_git_add_res EQUAL 0)
  message(FATAL_ERROR "run_module_fetch_tests.cmake: git add failed in fixture dir (result ${_gsmf_git_add_res})")
endif()

execute_process(COMMAND "${_gsmf_git}"
    "--git-dir=${GSMF_FIXTURE_DIR}/.git" "--work-tree=${GSMF_FIXTURE_DIR}"
    -c user.name=gismo-test -c user.email=test@invalid -c commit.gpgsign=false
    commit -q -m fixture
  RESULT_VARIABLE _gsmf_git_commit_res
  OUTPUT_QUIET ERROR_QUIET)
if(NOT _gsmf_git_commit_res EQUAL 0)
  message(FATAL_ERROR "run_module_fetch_tests.cmake: git commit failed in fixture dir (result ${_gsmf_git_commit_res})")
endif()

set(GSMF_FIXTURE_URL "file://${GSMF_FIXTURE_DIR}")
## Common -D flags every case's sub-configure needs: GISMO_REPO=git skips
## the svn/zip-download branches of gismo_fetch_module_source() entirely
## (it only takes the "x${GISMO_REPO}" STREQUAL "xgit" branch),
## GISMO_FETCH_PROT=https pre-empts its own git-remote-based auto-detection,
## and <module>_url pins the clone source to the local fixture repository.
set(GSMF_COMMON_ARGS
  -DGISMO_CMAKE_DIR=${GISMO_SOURCE_DIR}/cmake
  -DGISMO_REPO=git
  -DGISMO_FETCH_PROT=https
  -DgsFetchFixture_url=${GSMF_FIXTURE_URL})

######################################################################
## case: absent_module
######################################################################
set(_case "absent_module")
gsmf_case_begin(${_case})
set(_case_dir "${SCRATCH_DIR}/${_case}")
file(REMOVE_RECURSE "${_case_dir}")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-src")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-bin")
gsmf_configure(${_case} "${GSMF_PROJECT_DIR}" "${_case_dir}/bin" _r _o
  -DCASE=${_case}
  -DSANDBOX_SRC_DIR=${_case_dir}/sandbox-src
  -DSANDBOX_BIN_DIR=${_case_dir}/sandbox-bin
  ${GSMF_COMMON_ARGS})
if(NOT _r EQUAL 0)
  message(STATUS "FAIL ${_case}/configure: expected '0' got '${_r}' (${_o})")
  _gsmf_inc(gsmf_fail_count)
else()
  gsmf_read_report("${_case_dir}/bin/module_fetch_report.txt" R)
  expect_bool(${_case} exists_before FALSE "${R_EXISTS_BEFORE}")
  expect_bool(${_case} cloned TRUE "${R_CLONED}")
  get_filename_component(_exp_inc "${_case_dir}/sandbox-src/optional/gsFetchFixture/include" REALPATH)
  get_filename_component(_exp_sentinel "${_case_dir}/sandbox-src/sentinel_include" REALPATH)
  expect_list_contains(${_case} include_dirs_final_fixture "${_exp_inc}" "${R_INCLUDE_DIRS_FINAL}")
  expect_list_contains(${_case} include_dirs_final_sentinel "${_exp_sentinel}" "${R_INCLUDE_DIRS_FINAL}")
  expect_str(${_case} deps_included "1" "${R_DEPS_INCLUDED}")
  expect_count_matches(${_case} deps_message "Processing dependencies for optional module: gsFetchFixture" 1 "${_o}")
endif()
gsmf_check_containment(${_case})
gsmf_case_end(${_case})

######################################################################
## case: legacy_order - the falsifier. Replays the OLD top-level order
## (dependencies hook, then gismo_fetch_module()) against the same fixture;
## the hook must find nothing, because the module has not been fetched yet.
######################################################################
set(_case "legacy_order")
gsmf_case_begin(${_case})
set(_case_dir "${SCRATCH_DIR}/${_case}")
file(REMOVE_RECURSE "${_case_dir}")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-src")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-bin")
gsmf_configure(${_case} "${GSMF_PROJECT_DIR}" "${_case_dir}/bin" _r _o
  -DCASE=${_case}
  -DSANDBOX_SRC_DIR=${_case_dir}/sandbox-src
  -DSANDBOX_BIN_DIR=${_case_dir}/sandbox-bin
  ${GSMF_COMMON_ARGS})
if(NOT _r EQUAL 0)
  message(STATUS "FAIL ${_case}/configure: expected '0' got '${_r}' (${_o})")
  _gsmf_inc(gsmf_fail_count)
else()
  gsmf_read_report("${_case_dir}/bin/module_fetch_report.txt" R)
  expect_bool(${_case} exists_before FALSE "${R_EXISTS_BEFORE}")
  get_filename_component(_exp_inc "${_case_dir}/sandbox-src/optional/gsFetchFixture/include" REALPATH)
  get_filename_component(_exp_sentinel "${_case_dir}/sandbox-src/sentinel_include" REALPATH)
  expect_list_contains(${_case} include_dirs_after_hook_sentinel "${_exp_sentinel}" "${R_INCLUDE_DIRS_AFTER_HOOK}")
  expect_list_not_contains(${_case} include_dirs_after_hook_fixture "${_exp_inc}" "${R_INCLUDE_DIRS_AFTER_HOOK}")
  expect_bool(${_case} cloned TRUE "${R_CLONED}")
  expect_str(${_case} deps_included "0" "${R_DEPS_INCLUDED}")
endif()
gsmf_check_containment(${_case})
gsmf_case_end(${_case})

######################################################################
## case: idempotent - two consecutive gismo_prepare_optional_modules() calls
## (plus a further gismo_fetch_module()) must neither re-clone nor re-include.
######################################################################
set(_case "idempotent")
gsmf_case_begin(${_case})
set(_case_dir "${SCRATCH_DIR}/${_case}")
file(REMOVE_RECURSE "${_case_dir}")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-src")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-bin")
gsmf_configure(${_case} "${GSMF_PROJECT_DIR}" "${_case_dir}/bin" _r _o
  -DCASE=${_case}
  -DSANDBOX_SRC_DIR=${_case_dir}/sandbox-src
  -DSANDBOX_BIN_DIR=${_case_dir}/sandbox-bin
  -DGISMO_SUBMODULES_HEAD=ON
  ${GSMF_COMMON_ARGS})
if(NOT _r EQUAL 0)
  message(STATUS "FAIL ${_case}/configure: expected '0' got '${_r}' (${_o})")
  _gsmf_inc(gsmf_fail_count)
else()
  gsmf_read_report("${_case_dir}/bin/module_fetch_report.txt" R)
  expect_bool(${_case} cloned TRUE "${R_CLONED}")
  expect_bool(${_case} marker_survives TRUE "${R_MARKER_SURVIVES}")
  expect_str(${_case} deps_included "1" "${R_DEPS_INCLUDED}")
  get_filename_component(_exp_inc "${_case_dir}/sandbox-src/optional/gsFetchFixture/include" REALPATH)
  expect_list_contains(${_case} include_dirs_final_fixture "${_exp_inc}" "${R_INCLUDE_DIRS_FINAL}")
  expect_count_matches(${_case} cloning_message "Cloning into gsFetchFixture" 1 "${_o}")
  expect_count_matches(${_case} deps_message "Processing dependencies for optional module: gsFetchFixture" 1 "${_o}")
endif()
gsmf_check_containment(${_case})
gsmf_case_end(${_case})

######################################################################
## summary - never FATAL_ERROR mid-run; every case above has already run
## regardless of earlier failures.
######################################################################
get_property(_final_case_fail GLOBAL PROPERTY gsmf_case_fail_count)
get_property(_final_pass GLOBAL PROPERTY gsmf_pass_count)
get_property(_final_skip GLOBAL PROPERTY gsmf_skip_count)
message(STATUS "RESULT: ${_final_pass} passed, ${_final_case_fail} failed, ${_final_skip} skipped")
if(_final_case_fail GREATER 0)
  message(FATAL_ERROR "run_module_fetch_tests.cmake: ${_final_case_fail} case(s) failed")
endif()
