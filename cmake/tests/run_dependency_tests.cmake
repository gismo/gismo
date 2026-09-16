######################################################################
## run_dependency_tests.cmake
##
## Standalone driver (invoked as `cmake -P run_dependency_tests.cmake`, wired
## into ctest as the single test `cmake_dependency_test`) for
## gismo_add_dependency() (cmake/gsFetch.cmake). Configures a battery of tiny
## fixture projects under a scratch directory below the CMake *binary* dir
## and asserts specific, individually-checkable observables from the
## find-or-fetch contract - target existence/type/include dirs, <Name>_FOUND/
## _VENDORED, GISMO_INCLUDE_DIRS/gismo_LINKER/gismo_EXTENSIONS before-vs-after,
## <Name>_OBJECTS, the compiled OBJECT library's type/visibility/compile
## options, exported dynamic symbols (via `nm -D` on a shared library built
## from those objects), and (for the NEVER cases) the exact failure text -
## never just "the command ran".
##
## Every case writes only under ${SCRATCH_DIR}; the containment assertions at
## the end of each case additionally check that nothing was written into the
## real source tree.
######################################################################

cmake_minimum_required(VERSION 3.1...3.10)

foreach(_gsdep_required GISMO_SOURCE_DIR SCRATCH_DIR TEST_GENERATOR)
  if(NOT DEFINED ${_gsdep_required})
    message(FATAL_ERROR "run_dependency_tests.cmake: -D${_gsdep_required}=... is required")
  endif()
endforeach()
if(NOT DEFINED GSDEP_TEST_NETWORK)
  set(GSDEP_TEST_NETWORK "AUTO")
endif()
## `cmake --build --config` needs a concrete configuration name even on a
## single-config generator that leaves TEST_BUILD_TYPE empty (a valid,
## common invocation) - default to Release so every `cmake --build` call
## below can pass `--config` unconditionally, which multi-config generators
## (Visual Studio, Ninja Multi-Config) require and single-config generators
## silently ignore.
if(NOT TEST_BUILD_TYPE)
  set(_gsdep_build_config "Release")
else()
  set(_gsdep_build_config "${TEST_BUILD_TYPE}")
endif()

set(GSDEP_TESTS_DIR    "${GISMO_SOURCE_DIR}/cmake/tests")
set(GSDEP_FIXTURES_DIR "${GSDEP_TESTS_DIR}/fixtures")
set(GSDEP_PROJECT_DIR  "${GSDEP_TESTS_DIR}/dependency_project")
set(GSDEP_CONSUMER_DIR "${GSDEP_TESTS_DIR}/consumer")

## Single source of truth for the subdirectory gismo_fetch_directory() fetches
## into (gsFetch.cmake's gismo_add_dependency(): `DESTINATION external`).
set(VENDOR_DEST_SUBDIR "external")

file(REMOVE_RECURSE "${SCRATCH_DIR}")
file(MAKE_DIRECTORY "${SCRATCH_DIR}")

set_property(GLOBAL PROPERTY gsdep_fail_count 0)
set_property(GLOBAL PROPERTY gsdep_pass_count 0)
set_property(GLOBAL PROPERTY gsdep_skip_count 0)
set_property(GLOBAL PROPERTY gsdep_case_start_fail 0)
## gsdep_fail_count counts failing *assertions* (an already-FAIL case can
## still print several); gsdep_case_fail_count counts failing *cases* and is
## what the final RESULT: line and the FATAL_ERROR gate use, so "<m> failed"
## always matches the number of "CASE <id> FAIL" lines above it.
set_property(GLOBAL PROPERTY gsdep_case_fail_count 0)

######################################################################
## expectation helpers - each prints exactly
##   FAIL <case>/<label>: expected '<expected>' got '<actual>'
## on mismatch and increments the shared fail counter (a GLOBAL PROPERTY,
## since these run inside function() scopes nested arbitrarily deep).
######################################################################

function(_gsdep_inc PROP)
  get_property(_v GLOBAL PROPERTY ${PROP})
  math(EXPR _v "${_v}+1")
  set_property(GLOBAL PROPERTY ${PROP} ${_v})
endfunction()

function(expect_str CASE_ID LABEL EXPECTED ACTUAL)
  if(NOT "${EXPECTED}" STREQUAL "${ACTUAL}")
    message(STATUS "FAIL ${CASE_ID}/${LABEL}: expected '${EXPECTED}' got '${ACTUAL}'")
    _gsdep_inc(gsdep_fail_count)
  endif()
endfunction()

function(expect_str_ne CASE_ID LABEL FORBIDDEN ACTUAL)
  # ACTUAL must differ from FORBIDDEN and must not be the "<unset>" sentinel
  # (an absent value can never satisfy "is not X" in any useful sense).
  if("${ACTUAL}" STREQUAL "${FORBIDDEN}" OR "${ACTUAL}" STREQUAL "<unset>")
    message(STATUS "FAIL ${CASE_ID}/${LABEL}: expected 'not ${FORBIDDEN} (and not <unset>)' got '${ACTUAL}'")
    _gsdep_inc(gsdep_fail_count)
  endif()
endfunction()

function(expect_match CASE_ID LABEL PATTERN ACTUAL)
  # CMake's own message(FATAL_ERROR ...)/message(SEND_ERROR ...) word-wraps
  # long text at a column width when it prints to the console, so a literal
  # phrase from gsFetch.cmake's error text can arrive here split across a
  # line break and re-indented. Collapse every whitespace run to one space
  # before matching so PATTERN only has to describe the words, not the
  # console's incidental wrapping.
  string(REGEX REPLACE "[ \t\r\n]+" " " _actual_norm "${ACTUAL}")
  if(NOT "${_actual_norm}" MATCHES "${PATTERN}")
    message(STATUS "FAIL ${CASE_ID}/${LABEL}: expected '${PATTERN}' got '${ACTUAL}'")
    _gsdep_inc(gsdep_fail_count)
  endif()
endfunction()

function(expect_bool CASE_ID LABEL EXPECTED_BOOL ACTUAL)
  # EXPECTED_BOOL is the literal 'TRUE' or 'FALSE'. ACTUAL is the raw report
  # value: TRUE/ON/1/YES are true, FALSE/OFF/0/NO/empty are false, and the
  # literal '<unset>' is always a FAIL - including when EXPECTED_BOOL is
  # FALSE, since a boolean gismo_add_dependency() never set is not "false",
  # it is unknown, and that distinction is exactly what '<unset>' exists to
  # preserve (an empty-vs-empty comparison would otherwise pass silently).
  if("${ACTUAL}" STREQUAL "<unset>")
    message(STATUS "FAIL ${CASE_ID}/${LABEL}: expected '${EXPECTED_BOOL}' got '<unset>'")
    _gsdep_inc(gsdep_fail_count)
    return()
  endif()
  if("${ACTUAL}" MATCHES "^(TRUE|ON|1|YES)$")
    set(_actual_bool "TRUE")
  elseif("${ACTUAL}" STREQUAL "" OR "${ACTUAL}" MATCHES "^(FALSE|OFF|0|NO)$")
    set(_actual_bool "FALSE")
  else()
    set(_actual_bool "UNRECOGNISED(${ACTUAL})")
  endif()
  if(NOT "${_actual_bool}" STREQUAL "${EXPECTED_BOOL}")
    message(STATUS "FAIL ${CASE_ID}/${LABEL}: expected '${EXPECTED_BOOL}' got '${ACTUAL}'")
    _gsdep_inc(gsdep_fail_count)
  endif()
endfunction()

function(expect_list_contains CASE_ID LABEL EXPECTED_DIR ACTUAL_LIST)
  # ACTUAL_LIST is a report value: '<unset>', or a '|'-joined list in which
  # one element may still carry a $<BUILD_INTERFACE:...> wrapper (a correct
  # implementation may leave one there; get_target_property() returns it
  # literally) - stripped before comparing against the realpath'd expectation.
  if("${ACTUAL_LIST}" STREQUAL "<unset>")
    message(STATUS "FAIL ${CASE_ID}/${LABEL}: expected to contain '${EXPECTED_DIR}' got '<unset>'")
    _gsdep_inc(gsdep_fail_count)
    return()
  endif()
  string(REPLACE "|" ";" _elems "${ACTUAL_LIST}")
  set(_found FALSE)
  foreach(_e ${_elems})
    set(_stripped "${_e}")
    if(_stripped MATCHES "^\\$<BUILD_INTERFACE:(.*)>$")
      set(_stripped "${CMAKE_MATCH_1}")
    endif()
    if("${_stripped}" STREQUAL "${EXPECTED_DIR}")
      set(_found TRUE)
    endif()
  endforeach()
  if(NOT _found)
    message(STATUS "FAIL ${CASE_ID}/${LABEL}: expected '${EXPECTED_DIR}' got '${ACTUAL_LIST}'")
    _gsdep_inc(gsdep_fail_count)
  endif()
endfunction()

function(expect_not_match CASE_ID LABEL PATTERN ACTUAL)
  # Same whitespace normalisation as expect_match - FAIL when ACTUAL DOES
  # match PATTERN, the inverse assertion needed for e.g. "the hidden arm's
  # dynamic symbol table must not name the dependency's own symbol".
  string(REGEX REPLACE "[ \t\r\n]+" " " _actual_norm "${ACTUAL}")
  if("${_actual_norm}" MATCHES "${PATTERN}")
    message(STATUS "FAIL ${CASE_ID}/${LABEL}: expected NOT '${PATTERN}' got '${ACTUAL}'")
    _gsdep_inc(gsdep_fail_count)
  endif()
endfunction()

## Counts how many of PIPE_LIST's '|'-joined elements are exactly ELEMENT.
## list(FILTER) is CMake 3.6+, above this suite's 3.1 floor, so a foreach is
## used instead of relying on it.
function(gsdep_count_element PIPE_LIST ELEMENT OUT_VAR)
  string(REPLACE "|" ";" _elems "${PIPE_LIST}")
  set(_n 0)
  foreach(_e ${_elems})
    if("${_e}" STREQUAL "${ELEMENT}")
      math(EXPR _n "${_n}+1")
    endif()
  endforeach()
  set(${OUT_VAR} "${_n}" PARENT_SCOPE)
endfunction()

## Builds TARGET in BIN_DIR without running anything - the half of
## gsdep_build_and_run() the symbol-visibility cases need, since their
## artifact is a shared library inspected with `nm`, not an executable run
## for its stdout.
function(gsdep_build_target CASE_ID BIN_DIR TARGET)
  execute_process(
    COMMAND ${CMAKE_COMMAND} --build . --target ${TARGET} --config ${_gsdep_build_config}
    WORKING_DIRECTORY "${BIN_DIR}"
    RESULT_VARIABLE _bres
    OUTPUT_VARIABLE _bout
    ERROR_VARIABLE _bout)
  if(NOT _bres EQUAL 0)
    message(STATUS "FAIL ${CASE_ID}/build_${TARGET}: expected '0' got '${_bres}' (${_bout})")
    _gsdep_inc(gsdep_fail_count)
  endif()
endfunction()

######################################################################
## case bookkeeping and report parsing
######################################################################

function(gsdep_case_begin CASE_ID)
  get_property(_f GLOBAL PROPERTY gsdep_fail_count)
  set_property(GLOBAL PROPERTY gsdep_case_start_fail ${_f})
endfunction()

function(gsdep_case_end CASE_ID)
  get_property(_f0 GLOBAL PROPERTY gsdep_case_start_fail)
  get_property(_f1 GLOBAL PROPERTY gsdep_fail_count)
  if(_f1 GREATER _f0)
    message(STATUS "CASE ${CASE_ID} FAIL")
    _gsdep_inc(gsdep_case_fail_count)
  else()
    message(STATUS "CASE ${CASE_ID} PASS")
    _gsdep_inc(gsdep_pass_count)
  endif()
endfunction()

function(gsdep_read_report FILE OUT_PREFIX)
  file(STRINGS "${FILE}" _lines)
  foreach(_line ${_lines})
    if(_line MATCHES "^([A-Z_]+)=(.*)$")
      set(${OUT_PREFIX}_${CMAKE_MATCH_1} "${CMAKE_MATCH_2}" PARENT_SCOPE)
    endif()
  endforeach()
endfunction()

function(gsdep_check_containment CASE_ID DEP_NAME)
  if(EXISTS "${GISMO_SOURCE_DIR}/${VENDOR_DEST_SUBDIR}/${DEP_NAME}")
    message(STATUS "FAIL ${CASE_ID}/containment_external: expected 'NOT EXISTS' got 'EXISTS (${GISMO_SOURCE_DIR}/${VENDOR_DEST_SUBDIR}/${DEP_NAME})'")
    _gsdep_inc(gsdep_fail_count)
  endif()
  if(EXISTS "${GISMO_SOURCE_DIR}/cmake/tests/external/${DEP_NAME}")
    message(STATUS "FAIL ${CASE_ID}/containment_testsexternal: expected 'NOT EXISTS' got 'EXISTS (${GISMO_SOURCE_DIR}/cmake/tests/external/${DEP_NAME})'")
    _gsdep_inc(gsdep_fail_count)
  endif()
endfunction()

## Configures CASE_SRC_DIR into CASE_BIN_DIR, forwarding the generator/
## compiler/build-type triple plus every extra -D flag given as ARGN.
## CASE_BIN_DIR is created but never removed here - callers that need a
## clean directory call file(REMOVE_RECURSE) themselves first (every case
## does, except found_multi_twice's deliberate second call into the same
## directory).
function(gsdep_configure CASE_ID CASE_SRC_DIR CASE_BIN_DIR OUT_RESULT_VAR OUT_OUTPUT_VAR)
  file(MAKE_DIRECTORY "${CASE_BIN_DIR}")
  execute_process(
    COMMAND ${CMAKE_COMMAND}
      -G "${TEST_GENERATOR}"
      -DCMAKE_MAKE_PROGRAM=${TEST_MAKE_PROGRAM}
      -DCMAKE_CXX_COMPILER=${TEST_CXX_COMPILER}
      -DCMAKE_BUILD_TYPE=${TEST_BUILD_TYPE}
      ${ARGN}
      "${CASE_SRC_DIR}"
    WORKING_DIRECTORY "${CASE_BIN_DIR}"
    RESULT_VARIABLE _result
    OUTPUT_VARIABLE _output
    ERROR_VARIABLE _output)
  set(${OUT_RESULT_VAR} "${_result}" PARENT_SCOPE)
  set(${OUT_OUTPUT_VAR} "${_output}" PARENT_SCOPE)
endfunction()

## Builds TARGET_NAME in BIN_DIR, runs the resulting executable and checks
## its exit code is 0 and its stdout matches EXPECTED_STDOUT_PATTERN.
## `--config` is passed unconditionally (single-config generators ignore it,
## multi-config generators require it); the executable is looked up both
## with and without CMAKE_EXECUTABLE_SUFFIX since that variable is not
## available in `cmake -P` script mode (no project()) to build the name from.
function(gsdep_build_and_run CASE_ID BIN_DIR TARGET_NAME EXPECTED_STDOUT_PATTERN)
  execute_process(
    COMMAND ${CMAKE_COMMAND} --build . --target ${TARGET_NAME} --config ${_gsdep_build_config}
    WORKING_DIRECTORY "${BIN_DIR}"
    RESULT_VARIABLE _bres
    OUTPUT_VARIABLE _bout
    ERROR_VARIABLE _bout)
  if(NOT _bres EQUAL 0)
    message(STATUS "FAIL ${CASE_ID}/consumer_build: expected '0' got '${_bres}' (${_bout})")
    _gsdep_inc(gsdep_fail_count)
    return()
  endif()
  set(_exe "${BIN_DIR}/bin/${TARGET_NAME}")
  set(_exe_win "${_exe}.exe")
  if(EXISTS "${_exe}")
    set(_exe_found "${_exe}")
  elseif(EXISTS "${_exe_win}")
    set(_exe_found "${_exe_win}")
  else()
    message(STATUS "FAIL ${CASE_ID}/consumer_exists: expected 'EXISTS' got 'MISSING (${_exe} or ${_exe_win})'")
    _gsdep_inc(gsdep_fail_count)
    return()
  endif()
  execute_process(
    COMMAND "${_exe_found}"
    RESULT_VARIABLE _rres
    OUTPUT_VARIABLE _rout
    ERROR_VARIABLE _rout)
  expect_str(${CASE_ID} consumer_exit "0" "${_rres}")
  expect_match(${CASE_ID} consumer_stdout "${EXPECTED_STDOUT_PATTERN}" "${_rout}")
endfunction()

######################################################################
## prelude: build every "already installed" fixture package once, into a
## shared scratch prefix, before any found-based case configures. The two
## genex fixture directories are created here (never in the source tree) and
## reused, unmodified, by every genex_* case.
######################################################################

set(_found_prefix "${SCRATCH_DIR}/found")
set(_genex_dir_a "${SCRATCH_DIR}/genex_dirs/dirA")
set(_genex_dir_b "${SCRATCH_DIR}/genex_dirs/dirB")
file(MAKE_DIRECTORY "${_genex_dir_a}")
file(MAKE_DIRECTORY "${_genex_dir_b}")
get_filename_component(_genex_dir_a_real "${_genex_dir_a}" REALPATH)
get_filename_component(_genex_dir_b_real "${_genex_dir_b}" REALPATH)

set(_prebuilt_bin "${SCRATCH_DIR}/prebuilt-bin")
file(REMOVE_RECURSE "${_prebuilt_bin}")
file(MAKE_DIRECTORY "${_prebuilt_bin}")
execute_process(
  COMMAND ${CMAKE_COMMAND}
    -G "${TEST_GENERATOR}"
    -DCMAKE_MAKE_PROGRAM=${TEST_MAKE_PROGRAM}
    -DCMAKE_CXX_COMPILER=${TEST_CXX_COMPILER}
    -DCMAKE_BUILD_TYPE=${TEST_BUILD_TYPE}
    -DPREFIX_DIR=${_found_prefix}
    -DFIXTURES_DIR=${GSDEP_FIXTURES_DIR}
    -DDIRA=${_genex_dir_a}
    -DDIRB=${_genex_dir_b}
    "${GSDEP_FIXTURES_DIR}/prebuilt"
  WORKING_DIRECTORY "${_prebuilt_bin}"
  RESULT_VARIABLE _pre_cfg_result
  OUTPUT_VARIABLE _pre_out
  ERROR_VARIABLE _pre_out)
if(NOT _pre_cfg_result EQUAL 0)
  message(FATAL_ERROR "run_dependency_tests.cmake: prelude configure failed (${_pre_cfg_result}):\n${_pre_out}")
endif()
execute_process(
  COMMAND ${CMAKE_COMMAND} --build . --config ${_gsdep_build_config}
  WORKING_DIRECTORY "${_prebuilt_bin}"
  RESULT_VARIABLE _pre_build_result
  OUTPUT_VARIABLE _pre_bout
  ERROR_VARIABLE _pre_bout)
if(NOT _pre_build_result EQUAL 0)
  message(FATAL_ERROR "run_dependency_tests.cmake: prelude build failed (${_pre_build_result}):\n${_pre_bout}")
endif()

######################################################################
## case: header_found
######################################################################
set(_case "header_found")
gsdep_case_begin(${_case})
set(_case_dir "${SCRATCH_DIR}/${_case}")
file(REMOVE_RECURSE "${_case_dir}")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-src")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-bin")
gsdep_configure(${_case} "${GSDEP_PROJECT_DIR}" "${_case_dir}/bin" _r _o
  -DSANDBOX_SRC_DIR=${_case_dir}/sandbox-src
  -DSANDBOX_BIN_DIR=${_case_dir}/sandbox-bin
  -DGISMO_CMAKE_DIR=${GISMO_SOURCE_DIR}/cmake
  -DGSDEP_CONSUMER_DIR=${GSDEP_CONSUMER_DIR}
  -DCASE=${_case}
  -DGISMO_DEPENDENCY_FETCH=AUTO
  -DFakeHdrDep_DIR=${_found_prefix}/FakeHdrDep/lib/cmake/FakeHdrDep)
if(NOT _r EQUAL 0)
  message(STATUS "FAIL ${_case}/configure: expected '0' got '${_r}' (${_o})")
  _gsdep_inc(gsdep_fail_count)
else()
  gsdep_read_report("${_case_dir}/bin/dep_report.txt" R)
  expect_str(${_case} target_exists "TRUE" "${R_TARGET_EXISTS}")
  expect_bool(${_case} found TRUE "${R_FOUND}")
  expect_bool(${_case} vendored FALSE "${R_VENDORED}")
  get_filename_component(_exp_inc "${_found_prefix}/FakeHdrDep/include" REALPATH)
  expect_list_contains(${_case} target_include "${_exp_inc}" "${R_TARGET_INCLUDE}")
  expect_list_contains(${_case} include_dirs_after "${_exp_inc}" "${R_INCLUDE_DIRS_AFTER}")
  expect_str(${_case} include_dirs_before "<unset>" "${R_INCLUDE_DIRS_BEFORE}")
  expect_str(${_case} linker_unchanged "${R_LINKER_BEFORE}" "${R_LINKER_AFTER}")
  expect_str(${_case} objects_empty "" "${R_OBJECTS}")
  expect_str(${_case} extensions_after "<unset>" "${R_EXTENSIONS_AFTER}")
  gsdep_build_and_run(${_case} "${_case_dir}/bin" use_header "OK header 42")
endif()
gsdep_check_containment(${_case} FakeHdrDep)
gsdep_case_end(${_case})

######################################################################
## case: header_vendored
######################################################################
set(_case "header_vendored")
gsdep_case_begin(${_case})
set(_case_dir "${SCRATCH_DIR}/${_case}")
file(REMOVE_RECURSE "${_case_dir}")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-bin")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-src/${VENDOR_DEST_SUBDIR}")
file(COPY "${GSDEP_FIXTURES_DIR}/header_only/" DESTINATION "${_case_dir}/sandbox-src/${VENDOR_DEST_SUBDIR}/VendoredHdrDep")
gsdep_configure(${_case} "${GSDEP_PROJECT_DIR}" "${_case_dir}/bin" _r _o
  -DSANDBOX_SRC_DIR=${_case_dir}/sandbox-src
  -DSANDBOX_BIN_DIR=${_case_dir}/sandbox-bin
  -DGISMO_CMAKE_DIR=${GISMO_SOURCE_DIR}/cmake
  -DGSDEP_CONSUMER_DIR=${GSDEP_CONSUMER_DIR}
  -DCASE=${_case}
  -DGISMO_DEPENDENCY_FETCH=AUTO)
if(NOT _r EQUAL 0)
  message(STATUS "FAIL ${_case}/configure: expected '0' got '${_r}' (${_o})")
  _gsdep_inc(gsdep_fail_count)
else()
  gsdep_read_report("${_case_dir}/bin/dep_report.txt" R)
  expect_str(${_case} target_exists "TRUE" "${R_TARGET_EXISTS}")
  expect_bool(${_case} found TRUE "${R_FOUND}")
  expect_bool(${_case} vendored TRUE "${R_VENDORED}")
  get_filename_component(_exp_inc "${_case_dir}/sandbox-src/${VENDOR_DEST_SUBDIR}/VendoredHdrDep/include" REALPATH)
  expect_list_contains(${_case} target_include "${_exp_inc}" "${R_TARGET_INCLUDE}")
  expect_list_contains(${_case} include_dirs_after "${_exp_inc}" "${R_INCLUDE_DIRS_AFTER}")
  expect_str(${_case} include_dirs_before "<unset>" "${R_INCLUDE_DIRS_BEFORE}")
  expect_str(${_case} linker_unchanged "${R_LINKER_BEFORE}" "${R_LINKER_AFTER}")
  expect_str(${_case} objects_empty "" "${R_OBJECTS}")
  expect_str(${_case} extensions_after "<unset>" "${R_EXTENSIONS_AFTER}")
  if(EXISTS "${_case_dir}/sandbox-src/${VENDOR_DEST_SUBDIR}/VendoredHdrDep/VENDOR_FIXTURE_MARKER")
    set(_marker "PRESENT")
  else()
    set(_marker "MISSING")
  endif()
  expect_str(${_case} vendor_marker_survives "PRESENT" "${_marker}")
  gsdep_build_and_run(${_case} "${_case_dir}/bin" use_header "OK header 42")
endif()
gsdep_check_containment(${_case} VendoredHdrDep)
gsdep_case_end(${_case})

######################################################################
## case: header_never
######################################################################
set(_case "header_never")
gsdep_case_begin(${_case})
set(_case_dir "${SCRATCH_DIR}/${_case}")
file(REMOVE_RECURSE "${_case_dir}")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-src")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-bin")
gsdep_configure(${_case} "${GSDEP_PROJECT_DIR}" "${_case_dir}/bin" _r _o
  -DSANDBOX_SRC_DIR=${_case_dir}/sandbox-src
  -DSANDBOX_BIN_DIR=${_case_dir}/sandbox-bin
  -DGISMO_CMAKE_DIR=${GISMO_SOURCE_DIR}/cmake
  -DGSDEP_CONSUMER_DIR=${GSDEP_CONSUMER_DIR}
  -DCASE=${_case}
  -DGISMO_DEPENDENCY_FETCH=NEVER)
if(_r EQUAL 0)
  message(STATUS "FAIL ${_case}/configure: expected 'non-zero' got '0'")
  _gsdep_inc(gsdep_fail_count)
else()
  expect_match(${_case} error_names_dependency "MissingHdrDep" "${_o}")
  expect_match(${_case} error_names_knob "(GISMO_DEPENDENCY_FETCH|NEVER)" "${_o}")
endif()
if(EXISTS "${_case_dir}/sandbox-src/${VENDOR_DEST_SUBDIR}/MissingHdrDep")
  message(STATUS "FAIL ${_case}/nothing_created: expected 'NOT EXISTS' got 'EXISTS'")
  _gsdep_inc(gsdep_fail_count)
endif()
gsdep_check_containment(${_case} MissingHdrDep)
gsdep_case_end(${_case})

######################################################################
## case: sources_found
######################################################################
set(_case "sources_found")
gsdep_case_begin(${_case})
set(_case_dir "${SCRATCH_DIR}/${_case}")
file(REMOVE_RECURSE "${_case_dir}")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-src")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-bin")
gsdep_configure(${_case} "${GSDEP_PROJECT_DIR}" "${_case_dir}/bin" _r _o
  -DSANDBOX_SRC_DIR=${_case_dir}/sandbox-src
  -DSANDBOX_BIN_DIR=${_case_dir}/sandbox-bin
  -DGISMO_CMAKE_DIR=${GISMO_SOURCE_DIR}/cmake
  -DGSDEP_CONSUMER_DIR=${GSDEP_CONSUMER_DIR}
  -DCASE=${_case}
  -DGISMO_DEPENDENCY_FETCH=AUTO
  -DFakeSrcDep_DIR=${_found_prefix}/FakeSrcDep/lib/cmake/FakeSrcDep)
if(NOT _r EQUAL 0)
  message(STATUS "FAIL ${_case}/configure: expected '0' got '${_r}' (${_o})")
  _gsdep_inc(gsdep_fail_count)
else()
  gsdep_read_report("${_case_dir}/bin/dep_report.txt" R)
  expect_str(${_case} target_exists "TRUE" "${R_TARGET_EXISTS}")
  expect_bool(${_case} found TRUE "${R_FOUND}")
  expect_bool(${_case} vendored FALSE "${R_VENDORED}")
  get_filename_component(_exp_inc "${_found_prefix}/FakeSrcDep/include" REALPATH)
  expect_list_contains(${_case} target_include "${_exp_inc}" "${R_TARGET_INCLUDE}")
  expect_list_contains(${_case} include_dirs_after "${_exp_inc}" "${R_INCLUDE_DIRS_AFTER}")
  expect_str(${_case} include_dirs_before "<unset>" "${R_INCLUDE_DIRS_BEFORE}")
  expect_str(${_case} linker_before "<unset>" "${R_LINKER_BEFORE}")
  expect_match(${_case} linker_after_has_entry "FakeSrcDep|gsdeptest_prebuilt" "${R_LINKER_AFTER}")
  expect_str(${_case} objects_empty "" "${R_OBJECTS}")
  expect_str(${_case} extensions_after "<unset>" "${R_EXTENSIONS_AFTER}")
  gsdep_build_and_run(${_case} "${_case_dir}/bin" use_sources "OK sources 4242")
endif()
gsdep_check_containment(${_case} FakeSrcDep)
gsdep_case_end(${_case})

######################################################################
## case: sources_vendored
######################################################################
set(_case "sources_vendored")
gsdep_case_begin(${_case})
set(_case_dir "${SCRATCH_DIR}/${_case}")
file(REMOVE_RECURSE "${_case_dir}")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-bin")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-src/${VENDOR_DEST_SUBDIR}")
file(COPY "${GSDEP_FIXTURES_DIR}/sources/" DESTINATION "${_case_dir}/sandbox-src/${VENDOR_DEST_SUBDIR}/VendoredSrcDep")
gsdep_configure(${_case} "${GSDEP_PROJECT_DIR}" "${_case_dir}/bin" _r _o
  -DSANDBOX_SRC_DIR=${_case_dir}/sandbox-src
  -DSANDBOX_BIN_DIR=${_case_dir}/sandbox-bin
  -DGISMO_CMAKE_DIR=${GISMO_SOURCE_DIR}/cmake
  -DGSDEP_CONSUMER_DIR=${GSDEP_CONSUMER_DIR}
  -DCASE=${_case}
  -DGISMO_DEPENDENCY_FETCH=AUTO)
if(NOT _r EQUAL 0)
  message(STATUS "FAIL ${_case}/configure: expected '0' got '${_r}' (${_o})")
  _gsdep_inc(gsdep_fail_count)
else()
  gsdep_read_report("${_case_dir}/bin/dep_report.txt" R)
  expect_str(${_case} target_exists "TRUE" "${R_TARGET_EXISTS}")
  expect_bool(${_case} found TRUE "${R_FOUND}")
  expect_bool(${_case} vendored TRUE "${R_VENDORED}")
  get_filename_component(_exp_inc "${_case_dir}/sandbox-src/${VENDOR_DEST_SUBDIR}/VendoredSrcDep/include" REALPATH)
  expect_list_contains(${_case} target_include "${_exp_inc}" "${R_TARGET_INCLUDE}")
  expect_list_contains(${_case} include_dirs_after "${_exp_inc}" "${R_INCLUDE_DIRS_AFTER}")
  expect_str(${_case} include_dirs_before "<unset>" "${R_INCLUDE_DIRS_BEFORE}")
  expect_str(${_case} target_type "INTERFACE_LIBRARY" "${R_TARGET_TYPE}")
  expect_str(${_case} real_target_type "OBJECT_LIBRARY" "${R_REAL_TARGET_TYPE}")
  expect_str(${_case} linker_before "<unset>" "${R_LINKER_BEFORE}")
  expect_str(${_case} linker_after "<unset>" "${R_LINKER_AFTER}")
  expect_str(${_case} objects "$<TARGET_OBJECTS:gismo_dep_VendoredSrcDep>" "${R_OBJECTS}")
  expect_str(${_case} extensions_before "<unset>" "${R_EXTENSIONS_BEFORE}")
  expect_str(${_case} extensions_after "$<TARGET_OBJECTS:gismo_dep_VendoredSrcDep>" "${R_EXTENSIONS_AFTER}")
  expect_match(${_case} compile_options_suppress_warnings "(^|\\|)(-w|/W0)(\\||$)" "${R_REAL_TARGET_COMPILE_OPTIONS}")
  if(EXISTS "${_case_dir}/sandbox-src/${VENDOR_DEST_SUBDIR}/VendoredSrcDep/VENDOR_FIXTURE_MARKER")
    set(_marker "PRESENT")
  else()
    set(_marker "MISSING")
  endif()
  expect_str(${_case} vendor_marker_survives "PRESENT" "${_marker}")
  # The consumer links no library directly (the INTERFACE target carries only
  # the include dir); it can only resolve gsdeptest_lib_answer() because
  # gismo_EXTENSIONS - which now holds this dependency's objects - is passed
  # as one of its sources (dependency_project/CMakeLists.txt add_executable call).
  gsdep_build_and_run(${_case} "${_case_dir}/bin" use_sources "OK sources 4242")
endif()
gsdep_check_containment(${_case} VendoredSrcDep)
gsdep_case_end(${_case})

######################################################################
## case: sources_vendored_twice - same scratch build dir configured TWICE (no
## file(REMOVE_RECURSE) between), following the found_multi_twice pattern for
## a vendored OBJECT library: the second configure must find gismo_EXTENSIONS
## already caching a pre-existing entry plus run 1's own
## $<TARGET_OBJECTS:...> entry, and dedup against the latter rather than
## re-appending or overwriting either.
######################################################################
set(_case "sources_vendored_twice")
gsdep_case_begin(${_case})
set(_case_dir "${SCRATCH_DIR}/${_case}")
file(REMOVE_RECURSE "${_case_dir}")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-bin")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-src/${VENDOR_DEST_SUBDIR}")
file(COPY "${GSDEP_FIXTURES_DIR}/sources/" DESTINATION "${_case_dir}/sandbox-src/${VENDOR_DEST_SUBDIR}/TwiceSrcDep")
set(_common_args
  -DSANDBOX_SRC_DIR=${_case_dir}/sandbox-src
  -DSANDBOX_BIN_DIR=${_case_dir}/sandbox-bin
  -DGISMO_CMAKE_DIR=${GISMO_SOURCE_DIR}/cmake
  -DGSDEP_CONSUMER_DIR=${GSDEP_CONSUMER_DIR}
  -DCASE=${_case}
  -DGISMO_DEPENDENCY_FETCH=AUTO)
gsdep_configure(${_case} "${GSDEP_PROJECT_DIR}" "${_case_dir}/bin" _r1 _o1 ${_common_args}
  -Dgismo_EXTENSIONS:INTERNAL=${GSDEP_CONSUMER_DIR}/extension_seed.cpp)
if(NOT _r1 EQUAL 0)
  message(STATUS "FAIL ${_case}/configure_run1: expected '0' got '${_r1}' (${_o1})")
  _gsdep_inc(gsdep_fail_count)
else()
  gsdep_read_report("${_case_dir}/bin/dep_report.txt" R1)
  expect_str(${_case} target_exists_run1 "TRUE" "${R1_TARGET_EXISTS}")
  # Second configure into the SAME build dir, no seed re-passed and no
  # REMOVE_RECURSE - gismo_EXTENSIONS persists as CACHE INTERNAL from run 1.
  gsdep_configure(${_case} "${GSDEP_PROJECT_DIR}" "${_case_dir}/bin" _r2 _o2 ${_common_args})
  if(NOT _r2 EQUAL 0)
    message(STATUS "FAIL ${_case}/configure_run2: expected '0' got '${_r2}' (${_o2})")
    _gsdep_inc(gsdep_fail_count)
  else()
    gsdep_read_report("${_case_dir}/bin/dep_report.txt" R2)
    expect_str(${_case} target_exists_run2 "TRUE" "${R2_TARGET_EXISTS}")
    # Proves run 2 re-entered the function with the run-1 list already in the
    # cache, so the dedup path was exercised and did not merely append once
    # to an empty list.
    expect_str(${_case} extensions_before_run2 "${R1_EXTENSIONS_AFTER}" "${R2_EXTENSIONS_BEFORE}")
    expect_str(${_case} extensions_after_stable "${R1_EXTENSIONS_AFTER}" "${R2_EXTENSIONS_AFTER}")
    gsdep_count_element("${R2_EXTENSIONS_AFTER}" "$<TARGET_OBJECTS:gismo_dep_TwiceSrcDep>" _n_objects)
    expect_str(${_case} extensions_object_count "1" "${_n_objects}")
    string(REPLACE "|" ";" _r2_ext_elems "${R2_EXTENSIONS_AFTER}")
    list(LENGTH _r2_ext_elems _n_total)
    expect_str(${_case} extensions_total_count "2" "${_n_total}")
    list(GET _r2_ext_elems 0 _r2_ext_first)
    expect_str(${_case} extensions_seed_first "${GSDEP_CONSUMER_DIR}/extension_seed.cpp" "${_r2_ext_first}")
    expect_str(${_case} linker_after "<unset>" "${R2_LINKER_AFTER}")
    expect_str(${_case} real_target_type "OBJECT_LIBRARY" "${R2_REAL_TARGET_TYPE}")
    expect_str(${_case} objects "$<TARGET_OBJECTS:gismo_dep_TwiceSrcDep>" "${R2_OBJECTS}")
    gsdep_build_and_run(${_case} "${_case_dir}/bin" use_sources "OK sources 4242")
  endif()
endif()
gsdep_check_containment(${_case} TwiceSrcDep)
gsdep_case_end(${_case})

######################################################################
## case: sources_export_symbols - EXPORT_SYMBOLS on a vendored MODE SOURCES
## dependency must keep its symbols visible in a shared library that links
## its objects, even under a hidden-by-default preset (the fixture sets
## CMAKE_CXX_VISIBILITY_PRESET hidden etc. itself, mirroring
## cmake/gsConfig.cmake:20-22, which this standalone project never includes).
## ELF dynamic symbol table only - ruled out on Windows/macOS, and needs a
## toolchain `nm`; `nm -D` rather than a "must fail to link" check, since the
## latter passes vacuously on any unrelated build error, whereas the anchor
## symbol (present in the same table in both arms) is a positive control.
######################################################################
set(_case "sources_export_symbols")
set(_case_skipped FALSE)
gsdep_case_begin(${_case})
if(CMAKE_HOST_WIN32 OR CMAKE_HOST_APPLE)
  message(STATUS "SKIP: ${_case} (ELF dynamic symbol table only)")
  message(STATUS "CASE ${_case} SKIP")
  _gsdep_inc(gsdep_skip_count)
  set(_case_skipped TRUE)
else()
  set(_case_dir "${SCRATCH_DIR}/${_case}")
  file(REMOVE_RECURSE "${_case_dir}")
  file(MAKE_DIRECTORY "${_case_dir}/sandbox-bin")
  file(MAKE_DIRECTORY "${_case_dir}/sandbox-src/${VENDOR_DEST_SUBDIR}")
  file(COPY "${GSDEP_FIXTURES_DIR}/sources/" DESTINATION "${_case_dir}/sandbox-src/${VENDOR_DEST_SUBDIR}/ExportSrcDep")
  gsdep_configure(${_case} "${GSDEP_PROJECT_DIR}" "${_case_dir}/bin" _r _o
    -DSANDBOX_SRC_DIR=${_case_dir}/sandbox-src
    -DSANDBOX_BIN_DIR=${_case_dir}/sandbox-bin
    -DGISMO_CMAKE_DIR=${GISMO_SOURCE_DIR}/cmake
    -DGSDEP_CONSUMER_DIR=${GSDEP_CONSUMER_DIR}
    -DCASE=${_case}
    -DGISMO_DEPENDENCY_FETCH=AUTO)
  if(NOT _r EQUAL 0)
    message(STATUS "FAIL ${_case}/configure: expected '0' got '${_r}' (${_o})")
    _gsdep_inc(gsdep_fail_count)
  else()
    gsdep_read_report("${_case_dir}/bin/dep_report.txt" R)
    if(R_NM STREQUAL "<unset>")
      message(STATUS "SKIP: ${_case} (no CMAKE_NM in the fixture's toolchain)")
      message(STATUS "CASE ${_case} SKIP")
      _gsdep_inc(gsdep_skip_count)
      set(_case_skipped TRUE)
    else()
      expect_str(${_case} real_target_type "OBJECT_LIBRARY" "${R_REAL_TARGET_TYPE}")
      gsdep_build_target(${_case} "${_case_dir}/bin" gsdeptest_shared)
      if(NOT EXISTS "${R_SHARED_LIB}")
        message(STATUS "FAIL ${_case}/shared_lib_exists: expected 'EXISTS' got 'MISSING (${R_SHARED_LIB})'")
        _gsdep_inc(gsdep_fail_count)
      else()
        execute_process(
          COMMAND "${R_NM}" -D --defined-only "${R_SHARED_LIB}"
          RESULT_VARIABLE _nm_res
          OUTPUT_VARIABLE _nm_out
          ERROR_VARIABLE _nm_out)
        expect_str(${_case} nm_exit "0" "${_nm_res}")
        expect_match(${_case} nm_sees_anchor "gsdeptest_shared_anchor" "${_nm_out}")
        expect_match(${_case} nm_exports_dep_symbol "gsdeptest_lib_answer" "${_nm_out}")
        expect_str(${_case} real_target_visibility "default|default|OFF" "${R_REAL_TARGET_VISIBILITY}")
      endif()
    endif()
  endif()
  gsdep_check_containment(${_case} ExportSrcDep)
endif()
if(NOT _case_skipped)
  gsdep_case_end(${_case})
endif()

######################################################################
## case: sources_hidden_symbols - the mirror of sources_export_symbols
## without EXPORT_SYMBOLS: the dependency's own symbol must NOT reach the
## shared library's dynamic table, while the anchor (a positive control
## exported unconditionally by the fixture) still does - without that
## control, an empty `nm -D` output would make this arm pass vacuously.
######################################################################
set(_case "sources_hidden_symbols")
set(_case_skipped FALSE)
gsdep_case_begin(${_case})
if(CMAKE_HOST_WIN32 OR CMAKE_HOST_APPLE)
  message(STATUS "SKIP: ${_case} (ELF dynamic symbol table only)")
  message(STATUS "CASE ${_case} SKIP")
  _gsdep_inc(gsdep_skip_count)
  set(_case_skipped TRUE)
else()
  set(_case_dir "${SCRATCH_DIR}/${_case}")
  file(REMOVE_RECURSE "${_case_dir}")
  file(MAKE_DIRECTORY "${_case_dir}/sandbox-bin")
  file(MAKE_DIRECTORY "${_case_dir}/sandbox-src/${VENDOR_DEST_SUBDIR}")
  file(COPY "${GSDEP_FIXTURES_DIR}/sources/" DESTINATION "${_case_dir}/sandbox-src/${VENDOR_DEST_SUBDIR}/HiddenSrcDep")
  gsdep_configure(${_case} "${GSDEP_PROJECT_DIR}" "${_case_dir}/bin" _r _o
    -DSANDBOX_SRC_DIR=${_case_dir}/sandbox-src
    -DSANDBOX_BIN_DIR=${_case_dir}/sandbox-bin
    -DGISMO_CMAKE_DIR=${GISMO_SOURCE_DIR}/cmake
    -DGSDEP_CONSUMER_DIR=${GSDEP_CONSUMER_DIR}
    -DCASE=${_case}
    -DGISMO_DEPENDENCY_FETCH=AUTO)
  if(NOT _r EQUAL 0)
    message(STATUS "FAIL ${_case}/configure: expected '0' got '${_r}' (${_o})")
    _gsdep_inc(gsdep_fail_count)
  else()
    gsdep_read_report("${_case_dir}/bin/dep_report.txt" R)
    if(R_NM STREQUAL "<unset>")
      message(STATUS "SKIP: ${_case} (no CMAKE_NM in the fixture's toolchain)")
      message(STATUS "CASE ${_case} SKIP")
      _gsdep_inc(gsdep_skip_count)
      set(_case_skipped TRUE)
    else()
      expect_str(${_case} real_target_type "OBJECT_LIBRARY" "${R_REAL_TARGET_TYPE}")
      gsdep_build_target(${_case} "${_case_dir}/bin" gsdeptest_shared)
      if(NOT EXISTS "${R_SHARED_LIB}")
        message(STATUS "FAIL ${_case}/shared_lib_exists: expected 'EXISTS' got 'MISSING (${R_SHARED_LIB})'")
        _gsdep_inc(gsdep_fail_count)
      else()
        execute_process(
          COMMAND "${R_NM}" -D --defined-only "${R_SHARED_LIB}"
          RESULT_VARIABLE _nm_res
          OUTPUT_VARIABLE _nm_out
          ERROR_VARIABLE _nm_out)
        expect_str(${_case} nm_exit "0" "${_nm_res}")
        expect_match(${_case} nm_sees_anchor "gsdeptest_shared_anchor" "${_nm_out}")
        expect_not_match(${_case} nm_hides_dep_symbol "gsdeptest_lib_answer" "${_nm_out}")
        # Only the CXX component is asserted: the C and inlines components of
        # REAL_TARGET_VISIBILITY depend on property initialisation for a
        # language (C) this fixture never enables.
        expect_match(${_case} real_target_visibility_cxx_hidden "^hidden\\|" "${R_REAL_TARGET_VISIBILITY}")
      endif()
    endif()
  endif()
  gsdep_check_containment(${_case} HiddenSrcDep)
endif()
if(NOT _case_skipped)
  gsdep_case_end(${_case})
endif()

######################################################################
## case: header_export_symbols - EXPORT_SYMBOLS combined with MODE
## HEADER_ONLY is a FATAL_ERROR raised in argument validation: a HEADER_ONLY
## dependency compiles no code of its own, so there is no target whose
## symbol visibility could be set. The vendored fixture is present (the same
## header_only/ copy header_vendored uses), so the only possible source of
## the failure is the validation itself, not a missing include directory.
######################################################################
set(_case "header_export_symbols")
gsdep_case_begin(${_case})
set(_case_dir "${SCRATCH_DIR}/${_case}")
file(REMOVE_RECURSE "${_case_dir}")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-bin")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-src/${VENDOR_DEST_SUBDIR}")
file(COPY "${GSDEP_FIXTURES_DIR}/header_only/" DESTINATION "${_case_dir}/sandbox-src/${VENDOR_DEST_SUBDIR}/ExportHdrDep")
gsdep_configure(${_case} "${GSDEP_PROJECT_DIR}" "${_case_dir}/bin" _r _o
  -DSANDBOX_SRC_DIR=${_case_dir}/sandbox-src
  -DSANDBOX_BIN_DIR=${_case_dir}/sandbox-bin
  -DGISMO_CMAKE_DIR=${GISMO_SOURCE_DIR}/cmake
  -DGSDEP_CONSUMER_DIR=${GSDEP_CONSUMER_DIR}
  -DCASE=${_case}
  -DGISMO_DEPENDENCY_FETCH=AUTO)
if(_r EQUAL 0)
  message(STATUS "FAIL ${_case}/configure: expected 'non-zero' got '0'")
  _gsdep_inc(gsdep_fail_count)
else()
  expect_match(${_case} error_names_dependency "ExportHdrDep" "${_o}")
  expect_match(${_case} error_names_keyword "EXPORT_SYMBOLS" "${_o}")
  expect_match(${_case} error_names_mode "HEADER_ONLY" "${_o}")
endif()
gsdep_check_containment(${_case} ExportHdrDep)
gsdep_case_end(${_case})

######################################################################
## case: sources_found_export_symbols - EXPORT_SYMBOLS on a *found* (not
## vendored) MODE SOURCES dependency is a documented no-op: whether
## find_package() succeeds is machine-dependent, and once it has, no
## gismo_dep_<Name> target exists for this keyword to configure.
######################################################################
set(_case "sources_found_export_symbols")
gsdep_case_begin(${_case})
set(_case_dir "${SCRATCH_DIR}/${_case}")
file(REMOVE_RECURSE "${_case_dir}")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-src")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-bin")
gsdep_configure(${_case} "${GSDEP_PROJECT_DIR}" "${_case_dir}/bin" _r _o
  -DSANDBOX_SRC_DIR=${_case_dir}/sandbox-src
  -DSANDBOX_BIN_DIR=${_case_dir}/sandbox-bin
  -DGISMO_CMAKE_DIR=${GISMO_SOURCE_DIR}/cmake
  -DGSDEP_CONSUMER_DIR=${GSDEP_CONSUMER_DIR}
  -DCASE=${_case}
  -DGISMO_DEPENDENCY_FETCH=AUTO
  -DFakeSrcDep_DIR=${_found_prefix}/FakeSrcDep/lib/cmake/FakeSrcDep)
if(NOT _r EQUAL 0)
  message(STATUS "FAIL ${_case}/configure: expected '0' got '${_r}' (${_o})")
  _gsdep_inc(gsdep_fail_count)
else()
  gsdep_read_report("${_case_dir}/bin/dep_report.txt" R)
  expect_str(${_case} target_exists "TRUE" "${R_TARGET_EXISTS}")
  expect_bool(${_case} found TRUE "${R_FOUND}")
  expect_bool(${_case} vendored FALSE "${R_VENDORED}")
  expect_str(${_case} objects_empty "" "${R_OBJECTS}")
  expect_str(${_case} extensions_after "<unset>" "${R_EXTENSIONS_AFTER}")
  expect_str(${_case} real_target_type "<unset>" "${R_REAL_TARGET_TYPE}")
  expect_match(${_case} linker_after_has_entry "gsdeptest_prebuilt" "${R_LINKER_AFTER}")
  gsdep_build_and_run(${_case} "${_case_dir}/bin" use_sources "OK sources 4242")
endif()
gsdep_check_containment(${_case} FakeSrcDep)
gsdep_case_end(${_case})

######################################################################
## case: sources_never
######################################################################
set(_case "sources_never")
gsdep_case_begin(${_case})
set(_case_dir "${SCRATCH_DIR}/${_case}")
file(REMOVE_RECURSE "${_case_dir}")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-src")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-bin")
gsdep_configure(${_case} "${GSDEP_PROJECT_DIR}" "${_case_dir}/bin" _r _o
  -DSANDBOX_SRC_DIR=${_case_dir}/sandbox-src
  -DSANDBOX_BIN_DIR=${_case_dir}/sandbox-bin
  -DGISMO_CMAKE_DIR=${GISMO_SOURCE_DIR}/cmake
  -DGSDEP_CONSUMER_DIR=${GSDEP_CONSUMER_DIR}
  -DCASE=${_case}
  -DGISMO_DEPENDENCY_FETCH=NEVER)
if(_r EQUAL 0)
  message(STATUS "FAIL ${_case}/configure: expected 'non-zero' got '0'")
  _gsdep_inc(gsdep_fail_count)
else()
  expect_match(${_case} error_names_dependency "MissingSrcDep" "${_o}")
  expect_match(${_case} error_names_knob "(GISMO_DEPENDENCY_FETCH|NEVER)" "${_o}")
endif()
if(EXISTS "${_case_dir}/sandbox-src/${VENDOR_DEST_SUBDIR}/MissingSrcDep")
  message(STATUS "FAIL ${_case}/nothing_created: expected 'NOT EXISTS' got 'EXISTS'")
  _gsdep_inc(gsdep_fail_count)
endif()
gsdep_check_containment(${_case} MissingSrcDep)
gsdep_case_end(${_case})

######################################################################
## case: found_multi_twice - same scratch build dir configured TWICE (no
## file(REMOVE_RECURSE) between), asserting the element count *and* contents
## of GISMO_INCLUDE_DIRS/gismo_LINKER are identical after the second
## configure - a found package reporting multi-element legacy
## *_INCLUDE_DIRS/*_LIBRARIES must not grow either global on re-configure.
######################################################################
set(_case "found_multi_twice")
gsdep_case_begin(${_case})
set(_case_dir "${SCRATCH_DIR}/${_case}")
file(REMOVE_RECURSE "${_case_dir}")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-src")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-bin")
set(_common_args
  -DSANDBOX_SRC_DIR=${_case_dir}/sandbox-src
  -DSANDBOX_BIN_DIR=${_case_dir}/sandbox-bin
  -DGISMO_CMAKE_DIR=${GISMO_SOURCE_DIR}/cmake
  -DGSDEP_CONSUMER_DIR=${GSDEP_CONSUMER_DIR}
  -DCASE=${_case}
  -DGISMO_DEPENDENCY_FETCH=AUTO
  -DMultiTwiceDep_DIR=${_found_prefix}/MultiTwiceDep/lib/cmake/MultiTwiceDep)
gsdep_configure(${_case} "${GSDEP_PROJECT_DIR}" "${_case_dir}/bin" _r1 _o1 ${_common_args})
if(NOT _r1 EQUAL 0)
  message(STATUS "FAIL ${_case}/configure_run1: expected '0' got '${_r1}' (${_o1})")
  _gsdep_inc(gsdep_fail_count)
else()
  gsdep_read_report("${_case_dir}/bin/dep_report.txt" R1)
  expect_str(${_case} target_exists_run1 "TRUE" "${R1_TARGET_EXISTS}")
  # Second configure into the SAME build dir - no REMOVE_RECURSE - so
  # GISMO_INCLUDE_DIRS/gismo_LINKER persist as CACHE INTERNAL from run 1.
  gsdep_configure(${_case} "${GSDEP_PROJECT_DIR}" "${_case_dir}/bin" _r2 _o2 ${_common_args})
  if(NOT _r2 EQUAL 0)
    message(STATUS "FAIL ${_case}/configure_run2: expected '0' got '${_r2}' (${_o2})")
    _gsdep_inc(gsdep_fail_count)
  else()
    gsdep_read_report("${_case_dir}/bin/dep_report.txt" R2)
    expect_str(${_case} target_exists_run2 "TRUE" "${R2_TARGET_EXISTS}")
    expect_bool(${_case} found_run1 TRUE "${R1_FOUND}")
    expect_bool(${_case} vendored_run1 FALSE "${R1_VENDORED}")
    expect_str(${_case} include_dirs_stable "${R1_INCLUDE_DIRS_AFTER}" "${R2_INCLUDE_DIRS_AFTER}")
    expect_str(${_case} linker_stable "${R1_LINKER_AFTER}" "${R2_LINKER_AFTER}")
    # Guard against the "stable" comparisons above passing vacuously on two
    # empty strings: assert both dirs are genuinely present after run 2.
    get_filename_component(_exp_multi_a "${_found_prefix}/MultiTwiceDep/include_a" REALPATH)
    get_filename_component(_exp_multi_b "${_found_prefix}/MultiTwiceDep/include_b" REALPATH)
    expect_list_contains(${_case} include_dirs_after_has_a "${_exp_multi_a}" "${R2_INCLUDE_DIRS_AFTER}")
    expect_list_contains(${_case} include_dirs_after_has_b "${_exp_multi_b}" "${R2_INCLUDE_DIRS_AFTER}")
    expect_match(${_case} linker_after_has_two_entries "liba\\.fake\\|.*libb\\.fake|libb\\.fake\\|.*liba\\.fake" "${R2_LINKER_AFTER}")
  endif()
endif()
gsdep_check_containment(${_case} MultiTwiceDep)
gsdep_case_end(${_case})

######################################################################
## case: found_no_include - config-defined target with no
## INTERFACE_INCLUDE_DIRECTORIES property at all, so the found branch's
## fallback include-dir candidates are all empty; resolution must fail
## rather than silently caching a package as found with nothing resolved.
######################################################################
set(_case "found_no_include")
gsdep_case_begin(${_case})
set(_case_dir "${SCRATCH_DIR}/${_case}")
file(REMOVE_RECURSE "${_case_dir}")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-src")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-bin")
gsdep_configure(${_case} "${GSDEP_PROJECT_DIR}" "${_case_dir}/bin" _r _o
  -DSANDBOX_SRC_DIR=${_case_dir}/sandbox-src
  -DSANDBOX_BIN_DIR=${_case_dir}/sandbox-bin
  -DGISMO_CMAKE_DIR=${GISMO_SOURCE_DIR}/cmake
  -DGSDEP_CONSUMER_DIR=${GSDEP_CONSUMER_DIR}
  -DCASE=${_case}
  -DGISMO_DEPENDENCY_FETCH=AUTO
  -DNoIncludeDep_DIR=${_found_prefix}/NoIncludeDep/lib/cmake/NoIncludeDep)
if(_r EQUAL 0)
  message(STATUS "FAIL ${_case}/configure: expected 'non-zero' got '0'")
  _gsdep_inc(gsdep_fail_count)
else()
  expect_match(${_case} error_message "no include directory could be determined" "${_o}")
  expect_match(${_case} error_names_dependency "NoIncludeDep" "${_o}")
endif()
gsdep_check_containment(${_case} NoIncludeDep)
gsdep_case_end(${_case})

######################################################################
## case: genex_build_multi - one $<BUILD_INTERFACE:...> genex wrapping two
## paths; CMake list-splits it into two raw property elements before
## gismo_add_dependency() ever sees them, so both must be rejoined and
## resolved back into two entries, in order, on the target and in
## GISMO_INCLUDE_DIRS.
######################################################################
set(_case "genex_build_multi")
gsdep_case_begin(${_case})
set(_case_dir "${SCRATCH_DIR}/${_case}")
file(REMOVE_RECURSE "${_case_dir}")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-src")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-bin")
gsdep_configure(${_case} "${GSDEP_PROJECT_DIR}" "${_case_dir}/bin" _r _o
  -DSANDBOX_SRC_DIR=${_case_dir}/sandbox-src
  -DSANDBOX_BIN_DIR=${_case_dir}/sandbox-bin
  -DGISMO_CMAKE_DIR=${GISMO_SOURCE_DIR}/cmake
  -DGSDEP_CONSUMER_DIR=${GSDEP_CONSUMER_DIR}
  -DCASE=${_case}
  -DGISMO_DEPENDENCY_FETCH=AUTO
  -DGenexBuildMultiDep_DIR=${_found_prefix}/GenexBuildMultiDep/lib/cmake/GenexBuildMultiDep)
if(NOT _r EQUAL 0)
  message(STATUS "FAIL ${_case}/configure: expected '0' got '${_r}' (${_o})")
  _gsdep_inc(gsdep_fail_count)
else()
  gsdep_read_report("${_case_dir}/bin/dep_report.txt" R)
  expect_str(${_case} target_exists "TRUE" "${R_TARGET_EXISTS}")
  expect_bool(${_case} found TRUE "${R_FOUND}")
  expect_bool(${_case} vendored FALSE "${R_VENDORED}")
  expect_str(${_case} target_include "${_genex_dir_a_real}|${_genex_dir_b_real}" "${R_TARGET_INCLUDE}")
  expect_str(${_case} include_dirs_after "${_genex_dir_a_real}|${_genex_dir_b_real}" "${R_INCLUDE_DIRS_AFTER}")
endif()
gsdep_check_containment(${_case} GenexBuildMultiDep)
gsdep_case_end(${_case})

######################################################################
## case: genex_concat - $<BUILD_INTERFACE:A>$<INSTALL_INTERFACE:include>
## concatenated in one property element with no separator, as a single
## target_include_directories() call with both arguments produces; only the
## BUILD_INTERFACE payload must survive, never a garbage concatenation of
## the two.
######################################################################
set(_case "genex_concat")
gsdep_case_begin(${_case})
set(_case_dir "${SCRATCH_DIR}/${_case}")
file(REMOVE_RECURSE "${_case_dir}")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-src")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-bin")
gsdep_configure(${_case} "${GSDEP_PROJECT_DIR}" "${_case_dir}/bin" _r _o
  -DSANDBOX_SRC_DIR=${_case_dir}/sandbox-src
  -DSANDBOX_BIN_DIR=${_case_dir}/sandbox-bin
  -DGISMO_CMAKE_DIR=${GISMO_SOURCE_DIR}/cmake
  -DGSDEP_CONSUMER_DIR=${GSDEP_CONSUMER_DIR}
  -DCASE=${_case}
  -DGISMO_DEPENDENCY_FETCH=AUTO
  -DGenexConcatDep_DIR=${_found_prefix}/GenexConcatDep/lib/cmake/GenexConcatDep)
if(NOT _r EQUAL 0)
  message(STATUS "FAIL ${_case}/configure: expected '0' got '${_r}' (${_o})")
  _gsdep_inc(gsdep_fail_count)
else()
  gsdep_read_report("${_case_dir}/bin/dep_report.txt" R)
  expect_str(${_case} target_exists "TRUE" "${R_TARGET_EXISTS}")
  expect_bool(${_case} found TRUE "${R_FOUND}")
  expect_bool(${_case} vendored FALSE "${R_VENDORED}")
  expect_str(${_case} target_include "${_genex_dir_a_real}" "${R_TARGET_INCLUDE}")
  expect_str(${_case} include_dirs_after "${_genex_dir_a_real}" "${R_INCLUDE_DIRS_AFTER}")
endif()
gsdep_check_containment(${_case} GenexConcatDep)
gsdep_case_end(${_case})

######################################################################
## case: genex_install_only - only $<INSTALL_INTERFACE:include>; nothing
## survives normalisation, so resolution must fail.
######################################################################
set(_case "genex_install_only")
gsdep_case_begin(${_case})
set(_case_dir "${SCRATCH_DIR}/${_case}")
file(REMOVE_RECURSE "${_case_dir}")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-src")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-bin")
gsdep_configure(${_case} "${GSDEP_PROJECT_DIR}" "${_case_dir}/bin" _r _o
  -DSANDBOX_SRC_DIR=${_case_dir}/sandbox-src
  -DSANDBOX_BIN_DIR=${_case_dir}/sandbox-bin
  -DGISMO_CMAKE_DIR=${GISMO_SOURCE_DIR}/cmake
  -DGSDEP_CONSUMER_DIR=${GSDEP_CONSUMER_DIR}
  -DCASE=${_case}
  -DGISMO_DEPENDENCY_FETCH=AUTO
  -DGenexInstallOnlyDep_DIR=${_found_prefix}/GenexInstallOnlyDep/lib/cmake/GenexInstallOnlyDep)
if(_r EQUAL 0)
  message(STATUS "FAIL ${_case}/configure: expected 'non-zero' got '0'")
  _gsdep_inc(gsdep_fail_count)
else()
  expect_match(${_case} error_message "no include directory could be determined" "${_o}")
  expect_match(${_case} error_names_dependency "GenexInstallOnlyDep" "${_o}")
endif()
gsdep_check_containment(${_case} GenexInstallOnlyDep)
gsdep_case_end(${_case})

######################################################################
## case: genex_unknown - an unrecognised $<...> element alongside a real,
## plain path: the unrecognised element must survive, exempt from EXISTS,
## while the plain sibling is still checked and the STATUS line naming it
## appears exactly once.
######################################################################
set(_case "genex_unknown")
gsdep_case_begin(${_case})
set(_case_dir "${SCRATCH_DIR}/${_case}")
file(REMOVE_RECURSE "${_case_dir}")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-src")
file(MAKE_DIRECTORY "${_case_dir}/sandbox-bin")
gsdep_configure(${_case} "${GSDEP_PROJECT_DIR}" "${_case_dir}/bin" _r _o
  -DSANDBOX_SRC_DIR=${_case_dir}/sandbox-src
  -DSANDBOX_BIN_DIR=${_case_dir}/sandbox-bin
  -DGISMO_CMAKE_DIR=${GISMO_SOURCE_DIR}/cmake
  -DGSDEP_CONSUMER_DIR=${GSDEP_CONSUMER_DIR}
  -DCASE=${_case}
  -DGISMO_DEPENDENCY_FETCH=AUTO
  -DGenexUnknownDep_DIR=${_found_prefix}/GenexUnknownDep/lib/cmake/GenexUnknownDep)
if(NOT _r EQUAL 0)
  message(STATUS "FAIL ${_case}/configure: expected '0' got '${_r}' (${_o})")
  _gsdep_inc(gsdep_fail_count)
else()
  gsdep_read_report("${_case_dir}/bin/dep_report.txt" R)
  expect_str(${_case} target_exists "TRUE" "${R_TARGET_EXISTS}")
  expect_bool(${_case} found TRUE "${R_FOUND}")
  expect_bool(${_case} vendored FALSE "${R_VENDORED}")
  set(_exp_unknown_token "$<$<CONFIG:Debug>:${_genex_dir_b}>")
  expect_str(${_case} target_include "${_exp_unknown_token}|${_genex_dir_a_real}" "${R_TARGET_INCLUDE}")
  expect_str(${_case} include_dirs_after "${_exp_unknown_token}|${_genex_dir_a_real}" "${R_INCLUDE_DIRS_AFTER}")
  string(REGEX MATCHALL "is an unevaluated generator expression" _matches "${_o}")
  list(LENGTH _matches _n_matches)
  expect_str(${_case} status_line_count "1" "${_n_matches}")
endif()
gsdep_check_containment(${_case} GenexUnknownDep)
gsdep_case_end(${_case})

######################################################################
## case: header_network - the only case that touches the network. Probes
## once; skips visibly on no network (AUTO), fails loudly if the operator
## demanded network (GSDEP_TEST_NETWORK=ON) and it is not there, and never
## downgrades a post-probe clone failure to a skip.
######################################################################
set(_case "header_network")
gsdep_case_begin(${_case})
if(GSDEP_TEST_NETWORK STREQUAL "OFF")
  message(STATUS "SKIP: header_network (GSDEP_TEST_NETWORK=OFF)")
  message(STATUS "CASE header_network SKIP")
  _gsdep_inc(gsdep_skip_count)
else()
  file(DOWNLOAD "https://github.com/gismo/gsUnitTest/info/refs?service=git-upload-pack"
    "${SCRATCH_DIR}/probe.bin" TIMEOUT 15 INACTIVITY_TIMEOUT 15 STATUS _probe_status)
  list(GET _probe_status 0 _probe_code)
  list(GET _probe_status 1 _probe_msg)
  if(NOT _probe_code EQUAL 0)
    if(GSDEP_TEST_NETWORK STREQUAL "ON")
      message(STATUS "FAIL ${_case}/network_probe: expected '0' got '${_probe_code} (${_probe_msg})'")
      _gsdep_inc(gsdep_fail_count)
      gsdep_case_end(${_case})
    else()
      message(STATUS "SKIP: header_network (network probe failed: ${_probe_code} ${_probe_msg})")
      message(STATUS "CASE header_network SKIP")
      _gsdep_inc(gsdep_skip_count)
    endif()
  else()
    set(_case_dir "${SCRATCH_DIR}/${_case}")
    file(REMOVE_RECURSE "${_case_dir}")
    file(MAKE_DIRECTORY "${_case_dir}/sandbox-src")
    file(MAKE_DIRECTORY "${_case_dir}/sandbox-bin")
    gsdep_configure(${_case} "${GSDEP_PROJECT_DIR}" "${_case_dir}/bin" _r _o
      -DSANDBOX_SRC_DIR=${_case_dir}/sandbox-src
      -DSANDBOX_BIN_DIR=${_case_dir}/sandbox-bin
      -DGISMO_CMAKE_DIR=${GISMO_SOURCE_DIR}/cmake
      -DGSDEP_CONSUMER_DIR=${GSDEP_CONSUMER_DIR}
      -DCASE=${_case}
      -DGISMO_DEPENDENCY_FETCH=ALWAYS)
    if(NOT _r EQUAL 0)
      message(STATUS "FAIL ${_case}/configure: expected '0' got '${_r}' (${_o})")
      _gsdep_inc(gsdep_fail_count)
    else()
      gsdep_read_report("${_case_dir}/bin/dep_report.txt" R)
      expect_str(${_case} target_exists "TRUE" "${R_TARGET_EXISTS}")
      expect_bool(${_case} found TRUE "${R_FOUND}")
      expect_bool(${_case} vendored TRUE "${R_VENDORED}")
      set(_fetched_dir "${_case_dir}/sandbox-src/${VENDOR_DEST_SUBDIR}/NetHdrDep")
      if(EXISTS "${_fetched_dir}/UnitTestPP.h")
        set(_hdr_state "PRESENT")
      else()
        set(_hdr_state "MISSING")
      endif()
      expect_str(${_case} fetched_header "PRESENT" "${_hdr_state}")
      get_filename_component(_exp_inc "${_fetched_dir}" REALPATH)
      expect_list_contains(${_case} target_include "${_exp_inc}" "${R_TARGET_INCLUDE}")
      expect_list_contains(${_case} include_dirs_after "${_exp_inc}" "${R_INCLUDE_DIRS_AFTER}")
    endif()
    gsdep_check_containment(${_case} NetHdrDep)
    gsdep_case_end(${_case})
  endif()
endif()

######################################################################
## summary - never FATAL_ERROR mid-run; every case above has already run
## regardless of earlier failures.
######################################################################
get_property(_final_case_fail GLOBAL PROPERTY gsdep_case_fail_count)
get_property(_final_pass GLOBAL PROPERTY gsdep_pass_count)
get_property(_final_skip GLOBAL PROPERTY gsdep_skip_count)
message(STATUS "RESULT: ${_final_pass} passed, ${_final_case_fail} failed, ${_final_skip} skipped")
if(_final_case_fail GREATER 0)
  message(FATAL_ERROR "run_dependency_tests.cmake: ${_final_case_fail} case(s) failed")
endif()
