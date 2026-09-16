######################################################################
## gsFetch.cmake
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
######################################################################

include(CMakeParseArguments)

#note: latest CMake has FetchContent
function(gismo_fetch_directory)
  #use: gismo_fetch_directory(name GIT_REPOSITORY  ${git_repo} DESTINATION  optional)
  set(GF_NAME "${ARGV0}")
  set(oneValueArgs
    # Protect the following options
    DESTINATION
    SOURCE_DIR
    BINARY_DIR
    CONFIGURE_COMMAND
    BUILD_COMMAND
    INSTALL_COMMAND
    TEST_COMMAND )
  cmake_parse_arguments(GF "${GF_NAME}" "${oneValueArgs}" "" ${ARGN})
  #message( GF_UNPARSED_ARGUMENTS "= ${GF_UNPARSED_ARGUMENTS}")

  file(GLOB RESULT ${gismo_SOURCE_DIR}/${GF_DESTINATION}/${GF_NAME})
  list(LENGTH RESULT RESULT_LENGTH)
  if(NOT RESULT_LENGTH EQUAL 0)
    #message(STATUS "Fetch ${GF_DESTINATION} module ${GF_NAME} - found")
    return()
  endif()

  message(STATUS "Fetch ${GF_DESTINATION} module ${GF_NAME}")
  set(GF_SOURCE_DIR   "${gismo_SOURCE_DIR}/${GF_DESTINATION}/${GF_NAME}")
  set(GF_BINARY_DIR   "${gismo_BINARY_DIR}/${GF_DESTINATION}/${GF_NAME}_fetch")
  set(GF_DOWNLOAD_DIR "${gismo_BINARY_DIR}/${GF_DESTINATION}/${GF_NAME}_fetch")
  set(${GF_NAME}_SOURCE_DIR "${GF_SOURCE_DIR}" PARENT_SCOPE)
  set(${GF_NAME}_BINARY_DIR "${GF_BINARY_DIR}" PARENT_SCOPE)
  file(REMOVE "${GF_DOWNLOAD_DIR}/CMakeCache.txt")

  #  if(NOT EXISTS ${GF_DOWNLOAD_DIR}/CMakeLists.txt)
  file(WRITE ${GF_DOWNLOAD_DIR}/CMakeLists.txt "if(POLICY CMP0048)\ncmake_policy(SET CMP0048 NEW)\nendif()\nif(POLICY CMP0054)\ncmake_policy(SET CMP0054 NEW)\nendif()\nif(CMAKE_VERSION VERSION_LESS 3.19)\ncmake_minimum_required(VERSION 2.8.12)\nelse()\ncmake_minimum_required(VERSION 3.1...3.10)\nendif()\nproject(${GF_NAME}_fetch NONE)\ninclude(ExternalProject)\nExternalProject_Add(${GF_NAME}_fetch\n ${GF_UNPARSED_ARGUMENTS}\n SOURCE_DIR          \"${GF_SOURCE_DIR}\"\n BINARY_DIR          \"${GF_BINARY_DIR}\"\n CONFIGURE_COMMAND   \"\"\n BUILD_COMMAND       \"\"\n INSTALL_COMMAND     \"\"\n TEST_COMMAND        \"\"\n UPDATE_DISCONNECTED TRUE)\n")
  #  endif()

  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}"
    -D "CMAKE_MAKE_PROGRAM:FILE=${CMAKE_MAKE_PROGRAM}" .
    OUTPUT_QUIET
    RESULT_VARIABLE result
    WORKING_DIRECTORY "${GF_DOWNLOAD_DIR}" )
  if(result)
    message(SEND_ERROR "Configure step for ${GF_NAME} failed: ${result}")
  endif()

  # make sure that directory exists
  cmake_parse_arguments(GF "${GF_NAME}" "SVN_REPOSITORY" "" ${ARGN})
  if(DEFINED GF_SVN_REPOSITORY AND NOT EXISTS "${GF_SOURCE_DIR}/.svn")
    execute_process(COMMAND ${CMAKE_MAKE_PROGRAM} clean
      OUTPUT_QUIET
      WORKING_DIRECTORY "${GF_DOWNLOAD_DIR}" )
  endif()

  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    OUTPUT_QUIET
    RESULT_VARIABLE result
    WORKING_DIRECTORY "${GF_DOWNLOAD_DIR}" )
  if(result)
    message(SEND_ERROR "Build step for ${GF_NAME} failed: ${result}")
  endif()

  file(GLOB RESULT ${gismo_SOURCE_DIR}/${GF_DESTINATION}/${GF_NAME})
  list(LENGTH RESULT RESULT_LENGTH)
  if(RESULT_LENGTH EQUAL 0)
    message(SEND_ERROR "Fetch ${GF_DESTINATION} module ${GF_NAME} - not found")
  else()
    message(STATUS "Fetch ${GF_DESTINATION} module ${GF_NAME} - downloaded")
  endif()

endfunction()

function(gismo_add_extension SUBMODULE)
  if(TARGET ${SUBMODULE})
    return()
  endif()
  if(EXISTS "${gismo_SOURCE_DIR}/optional/${SUBMODULE}/CMakeLists.txt")
    add_subdirectory(${gismo_SOURCE_DIR}/optional/${SUBMODULE} ${gismo_BINARY_DIR}/optional/${SUBMODULE})
    if(EXISTS "${gismo_SOURCE_DIR}/optional/${SUBMODULE}/filedata")
      string(REGEX MATCH "optional/${SUBMODULE}/filedata" fmatch ${GISMO_SEARCH_PATHS})
      if(NOT fmatch)
        set(GISMO_SEARCH_PATHS "${GISMO_SEARCH_PATHS};${gismo_SOURCE_DIR}/optional/${SUBMODULE}/filedata/" CACHE INTERNAL "File search paths")
      endif()
    endif()
  else()
    message(WARNING "${SUBMODULE} does not contain CMakeLists.txt.")
  endif()

endfunction()

# called to fetch/download a submodule form git (working) and svn
# (ARGV0) SUBMODULE:  name of submodule
##
## gismo_fetch_module_source(SUBMODULE): fetch/update/restore a submodule's
## source tree only - never adds it as a build extension. Runs at most once
## per configure per module, latched on a GLOBAL property (not a cache
## variable, so a `cmake .` re-configure still re-runs the git update/restore
## - the same rationale as the gismo_add_dependency latch below).
function(gismo_fetch_module_source SUBMODULE)
  get_property(_gsmod_fetched GLOBAL PROPERTY GISMO_MODULE_FETCHED_${SUBMODULE} SET)
  if(_gsmod_fetched)
    return()
  endif()
  set_property(GLOBAL PROPERTY GISMO_MODULE_FETCHED_${SUBMODULE} TRUE)

  if(EXISTS "${gismo_SOURCE_DIR}/optional/${SUBMODULE}/CMakeLists.txt")
    #Update to current HEAD
    if(GISMO_SUBMODULES_HEAD AND EXISTS "${gismo_SOURCE_DIR}/optional/${SUBMODULE}/.git")
      message("Git fetch submodule ${SUBMODULE}")
      execute_process(COMMAND "${GIT_EXECUTABLE}" "fetch" "--depth" "1"
	ERROR_QUIET
	WORKING_DIRECTORY ${gismo_SOURCE_DIR}/optional/${SUBMODULE})
      execute_process(COMMAND "${GIT_EXECUTABLE}" "reset" "--hard"
	ERROR_QUIET
	WORKING_DIRECTORY ${gismo_SOURCE_DIR}/optional/${SUBMODULE})
      execute_process(COMMAND "${GIT_EXECUTABLE}" "clean" "-dfx"
	WORKING_DIRECTORY ${gismo_SOURCE_DIR}/optional/${SUBMODULE})
      execute_process(COMMAND "${GIT_EXECUTABLE}" "rebase"
	WORKING_DIRECTORY ${gismo_SOURCE_DIR}/optional/${SUBMODULE})
    endif()

    # Restore submodules to their proper hash
    # (!)Any local modifications will be LOST.)
    if(RESTORE_SUBMODULES AND NOT GISMO_SUBMODULES_HEAD AND EXISTS "${gismo_SOURCE_DIR}/optional/${SUBMODULE}/.git")
      message("Git restore ${SUBMODULE}")
      execute_process(COMMAND "${GIT_EXECUTABLE}" "restore"
	WORKING_DIRECTORY ${gismo_SOURCE_DIR}/optional/${SUBMODULE}
	ERROR_QUIET)
      if(${SUBMODULE}_HASH) # can be defined in submodules.txt
	execute_process(COMMAND "${GIT_EXECUTABLE}" "fetch" "--depth=1" origin ${${SUBMODULE}_HASH}
	  WORKING_DIRECTORY ${gismo_SOURCE_DIR}/optional/${SUBMODULE}
	  ERROR_QUIET
	  RESULT_VARIABLE gitfetch_res)
	if(gitfetch_res AND NOT gitfetch_res EQUAL 0)
	  message(FATAL_ERROR "Unable to fetch commit hash ${${SUBMODULE}_HASH} in ${SUBMODULE} (${git_repo})")
	endif()
	execute_process(COMMAND "${GIT_EXECUTABLE}" "checkout" ${${SUBMODULE}_HASH}
	  OUTPUT_QUIET
	  ERROR_QUIET
	  WORKING_DIRECTORY ${gismo_SOURCE_DIR}/optional/${SUBMODULE})
	execute_process(COMMAND "${GIT_EXECUTABLE}" "clean" "-dfx"
	  WORKING_DIRECTORY ${gismo_SOURCE_DIR}/optional/${SUBMODULE})
      endif()
      message("Hash is now ${${SUBMODULE}_HASH}")
    endif()

    # HERE:
    # add target unshallow_${SUBMODULE}

    return()
  endif()

  if (NOT DEFINED GISMO_FETCH_PROT)
    if(EXISTS "${gismo_SOURCE_DIR}/.git")
      execute_process(COMMAND "${GIT_EXECUTABLE}" "remote" "-v" OUTPUT_VARIABLE git_remote_res
	WORKING_DIRECTORY ${gismo_SOURCE_DIR})
      string(REGEX MATCH "git@github.com" fmatch "${git_remote_res}")
      if(fmatch)
	set(GISMO_FETCH_PROT "ssh" CACHE INTERNAL "")
      else()
	set(GISMO_FETCH_PROT "https" CACHE INTERNAL "")
      endif()
    endif()
    message(STATUS "Using ${GISMO_FETCH_PROT} git protocol")
  endif()

  get_repo_info(GISMO_REPO GISMO_REPO_REV) # or set manually
  #message("Fetch ${SUBMODULE} (repository: ${GISMO_REPO}, revision: ${GISMO_REPO_REV}, protocol: ${GISMO_FETCH_PROT}, username: ${GISMO_UNAME}, password: ${GISMO_PASS})")

  if("x${GISMO_REPO}" STREQUAL "xgit")
    if(NOT ${SUBMODULE}_url)
      if("x${GISMO_FETCH_PROT}" STREQUAL "xssh")
	set(${SUBMODULE}_url git@github.com:gismo/${SUBMODULE}.git)
      elseif("x${GISMO_FETCH_PROT}" STREQUAL "xhttps")
	set(${SUBMODULE}_url https://github.com/gismo/${SUBMODULE}.git)
      endif()
    endif()

    if(NOT EXISTS "${gismo_SOURCE_DIR}/optional/${SUBMODULE}/CMakeLists.txt")
      message(STATUS "Cloning into ${SUBMODULE}")
      find_package(Git REQUIRED)

      # Fetch SUBMODULE (note: git fetch --unshallow to get full clone)
      execute_process(COMMAND "${GIT_EXECUTABLE}" "clone" "--depth" "1" ${${SUBMODULE}_url}
	WORKING_DIRECTORY ${gismo_SOURCE_DIR}/optional
	ERROR_QUIET
	#OUTPUT_VARIABLE gitclone_out
	#ERROR_VARIABLE gitclone_err
	RESULT_VARIABLE gitclone_res)
      if(gitclone_res AND NOT gitclone_res EQUAL 0)
	message(FATAL_ERROR "Unable to clone module ${SUBMODULE} (${${SUBMODULE}_url}), set ${SUBMODULE}_url")
      endif()

      if(NOT GISMO_SUBMODULES_HEAD AND ${SUBMODULE}_HASH)# hash in submodules.txt
	execute_process(COMMAND "${GIT_EXECUTABLE}" "fetch" "--depth=1" origin ${${SUBMODULE}_HASH}
	  WORKING_DIRECTORY ${gismo_SOURCE_DIR}/optional/${SUBMODULE}
	  ERROR_QUIET
	  RESULT_VARIABLE gitfetch_res)
	if(gitfetch_res AND NOT gitfetch_res EQUAL 0)
	  message(FATAL_ERROR "Unable to fetch commit hash ${${SUBMODULE}_HASH} in ${SUBMODULE} (${${SUBMODULE}_url})")
	endif()
	execute_process(COMMAND "${GIT_EXECUTABLE}" "checkout" ${${SUBMODULE}_HASH}
	  OUTPUT_QUIET
	  ERROR_QUIET
	  WORKING_DIRECTORY ${gismo_SOURCE_DIR}/optional/${SUBMODULE})
	execute_process(COMMAND "${GIT_EXECUTABLE}" "clean" "-dfx"
	  WORKING_DIRECTORY ${gismo_SOURCE_DIR}/optional/${SUBMODULE})
	# note: specific commits are also reachable as, eg,
	# https://github.com/gismo/gismo/archive/<full-hash>.zip
      endif()
    endif()

  elseif("x${GISMO_REPO}" STREQUAL "xsvn")
    if("x${GISMO_FETCH_PROT}" STREQUAL "xssh")
      message(ERROR "GitHub does not support svn+ssh")
    endif()
    gismo_fetch_directory(${SUBMODULE}
      SVN_REPOSITORY https://github.com/gismo/${SUBMODULE}/trunk
      SVN_USERNAME ${GISMO_UNAME} # Username for Subversion checkout and update
      SVN_PASSWORD ${GISMO_PASS}  # Password for Subversion checkout and update
      SVN_TRUST_CERT 1            # Trust the Subversion server site certificate
      DESTINATION  optional )

  else()
    gismo_fetch_directory(${SUBMODULE}
      URL https://github.com/gismo/${SUBMODULE}/archive/master.zip
      DESTINATION  optional )
  endif()
endfunction()

## gismo_fetch_module(SUBMODULE): fetch the module's source
## (gismo_fetch_module_source, at most once per configure) and add it as a
## build extension (gismo_add_extension).
function(gismo_fetch_module SUBMODULE)
  gismo_fetch_module_source(${SUBMODULE})
  gismo_add_extension(${SUBMODULE})
endfunction()

## gismo_include_module_dependencies(SUBMODULE): include a module's
## dependencies.cmake exactly once, in the CALLER's scope. Must be a macro,
## not a function: a dependencies.cmake file sets plain variables (e.g.
## optional/gsAutoDiff/dependencies.cmake documents that it sets
## autodiff_FOUND and GISMO_INCLUDE_DIRS) that later configure logic reads -
## a function scope would swallow them. Latches only when the file is
## actually included, so a call made before the module's source exists does
## not block a later call once it does.
macro(gismo_include_module_dependencies SUBMODULE)
  get_property(_gsmod_deps_included GLOBAL PROPERTY GISMO_MODULE_DEPS_INCLUDED_${SUBMODULE} SET)
  if(NOT _gsmod_deps_included AND EXISTS "${gismo_SOURCE_DIR}/optional/${SUBMODULE}/dependencies.cmake")
    set_property(GLOBAL PROPERTY GISMO_MODULE_DEPS_INCLUDED_${SUBMODULE} TRUE)
    message(STATUS "Processing dependencies for optional module: ${SUBMODULE}")
    include("${gismo_SOURCE_DIR}/optional/${SUBMODULE}/dependencies.cmake")
  endif()
  unset(_gsmod_deps_included)
endmacro()

## gismo_prepare_optional_modules(): strip/dedup/sort GISMO_OPTIONAL, fetch
## the source of every listed module, then include every module's
## dependencies.cmake. Every source is fetched before any dependencies.cmake
## is included - needed when one module's dependencies.cmake reads or
## expects the source of another listed module, not just its own. Also a
## macro, so the fetch/include calls above and GISMO_OPTIONAL itself land in
## the caller's (top-level) scope.
macro(gismo_prepare_optional_modules)
  set(_gsmod_list "")
  foreach(_gsmod_m ${GISMO_OPTIONAL})
    string(STRIP "${_gsmod_m}" _gsmod_m)
    if(NOT _gsmod_m STREQUAL "")
      list(APPEND _gsmod_list "${_gsmod_m}")
    endif()
  endforeach()
  if(_gsmod_list)
    list(REMOVE_DUPLICATES _gsmod_list)
    list(SORT _gsmod_list)
  endif()
  set(GISMO_OPTIONAL ${_gsmod_list})
  foreach(_gsmod_m ${GISMO_OPTIONAL})
    gismo_fetch_module_source(${_gsmod_m})
  endforeach()
  foreach(_gsmod_m ${GISMO_OPTIONAL})
    gismo_include_module_dependencies(${_gsmod_m})
  endforeach()
  unset(_gsmod_m)
  unset(_gsmod_list)
endmacro()

######################################################################
## gismo_add_dependency: find-or-fetch resolution for an external dependency
##
## gismo_add_dependency(<Name>
##   [FIND_PACKAGE <pkg> [<version>] [COMPONENTS <c>...]]
##   [GIT_REPOSITORY <https url> [GIT_TAG <tag>] | URL <https url> [URL_HASH <hash>] | SOURCE_DIR <path>]
##   MODE HEADER_ONLY | SOURCES
##   [INCLUDE_SUBDIR <dir>]
##   [SOURCE_GLOBS <glob>...]
##   [EXPORT_SYMBOLS]
##   TARGET <Ns>::<Name>)
##
## Resolution order: if FIND_PACKAGE is given and the effective fetch mode is
## not ALWAYS, find_package(<pkg> ...) is tried first; on success an
## imported target is defined (or the target the found package's own config
## already created is reused) and nothing is fetched. Otherwise the
## dependency is fetched via gismo_fetch_directory() (unless SOURCE_DIR
## points at a local copy already) and vendored/compiled locally. If neither
## route can succeed the function stops the configure with FATAL_ERROR - it
## never leaves a partially-defined target behind.
##
## Two knobs steer this, both using the vocabulary AUTO | ALWAYS | NEVER:
##  - GISMO_DEPENDENCY_FETCH (cmake/gsOptions.cmake): the project-wide default.
##  - GISMO_<NAME>_FETCH (e.g. GISMO_SPECTRA_FETCH for gismo_add_dependency(Spectra ...)):
##    a per-dependency override, used when defined and non-empty.
## AUTO tries find_package() first and fetches on failure; ALWAYS always
## fetches/vendors, skipping find_package(); NEVER never fetches and fails
## the configure when find_package() does not succeed (or is not given).
## The legacy boolean GISMO_EIGEN_FETCH is a separate knob that this
## function does not consume. When SOURCE_DIR is given, nothing is fetched
## regardless of this knob, so the AUTO/ALWAYS/NEVER check below is skipped
## entirely on that branch - it only guards the GIT_REPOSITORY/URL path.
##
## Three CACHE INTERNAL outputs are set on success:
##  - <Name>_FOUND        TRUE (every failure path is FATAL_ERROR, so this
##                         function either resolves the dependency or stops
##                         the configure).
##  - <Name>_VENDORED     TRUE when fetched/compiled by G+Smo, FALSE when
##                         resolved via find_package() - GISMO_EIGEN_VENDORED
##                         (which gates header installation in
##                         cmake/gsInstall.cmake) is derived from
##                         EIGEN_INCLUDE_DIR's resolved path instead, so it no
##                         longer simply mirrors this flag.
##  - <Name>_OBJECTS      $<TARGET_OBJECTS:gismo_dep_<Name>> for a vendored
##                         MODE SOURCES dependency, empty string in every
##                         other case (found, or vendored HEADER_ONLY). It is
##                         the handle a per-module shared-library build would
##                         list among its own sources to receive the compiled
##                         dependency code, since such a build has no separate
##                         parameter for it (see the SOURCES bullet below).
##
## On the find branch (find_package() succeeded) TARGET is always an
## `add_library(... INTERFACE IMPORTED GLOBAL)`, regardless of MODE - MODE
## only steers what happens when nothing was found. A found dependency
## declared MODE HEADER_ONLY that nevertheless reports libraries (via
## <Pkg>_LIBRARY/_LIBRARIES, or read off TARGET as described below) still
## gets a gismo_LINKER entry: the found branch collects libraries
## unconditionally, MODE is not consulted there. The libraries themselves are
## resolved in two steps, the second only tried when the first finds nothing:
## first <Pkg>_LIBRARY/_LIBRARIES (any of the four case/plural spellings);
## then, when the package's own Config.cmake instead defines its imported
## target directly (add_library(...IMPORTED) + IMPORTED_LOCATION, the shape
## CMake itself recommends, and reports no such variable), TARGET's own link
## interface via `_gismo_dependency_link_paths()` (defined further down in
## this file) - TARGET's compiled-artifact location plus every absolute path,
## plain linker token, or nested imported target named in its
## INTERFACE_LINK_LIBRARIES, walked recursively with a cycle guard; a
## generator-expression entry cannot be evaluated at configure time and is
## skipped with a status message instead. On a multi-config generator
## (CMAKE_BUILD_TYPE empty) the artifact location this second step picks is
## necessarily the same for every configuration - accepted, since this
## function runs once at configure time, before a build type is selected.
## TARGET's own INTERFACE_LINK_LIBRARIES is set (via set_target_properties)
## from those libraries only when this function creates TARGET itself (the
## package's config did not already define one) - when the config already
## defined TARGET, that config's own INTERFACE_LINK_LIBRARIES is left
## untouched, whatever it is; the second resolution step only *reads* such a
## foreign target's properties, it never writes them.
## A found package's INTERFACE_INCLUDE_DIRECTORIES is normalised before use:
## $<BUILD_INTERFACE:X> unwraps to X, $<INSTALL_INTERFACE:X> is dropped
## (relative to the other project's install prefix, meaningless here), and
## any other $<...> is kept and exempted from the EXISTS check below, reported
## once via message(STATUS ...). If the package's config already defined
## TARGET, the normalised list is written back onto that (foreign) target's
## INTERFACE_INCLUDE_DIRECTORIES - deliberately: it is the only way TARGET can
## carry a usable build-tree path. The consequence: TARGET no longer carries
## whatever $<INSTALL_INTERFACE:...> entry the package's config gave it, so a
## later `install(EXPORT)` of a G+Smo target that links TARGET will export a
## target with no install-tree include path for that dependency.
## On the fetch/vendored branch the target shape follows MODE:
##  - HEADER_ONLY: TARGET is an `add_library(... INTERFACE IMPORTED GLOBAL)`
##    with INTERFACE_INCLUDE_DIRECTORIES set to the resolved include directory.
##  - SOURCES: the real compiled target is an OBJECT library named
##    `gismo_dep_<Name>` (never the caller-visible TARGET name, to avoid
##    collisions with module targets and with gismo_fetch_directory's own
##    naming), with a PRIVATE include dir, POSITION_INDEPENDENT_CODE ON, and
##    warnings suppressed (`-w` on GNU/Clang/AppleClang, `/W0` on MSVC) -
##    vendored code is not held to this project's own warning level.
##    SOURCE_GLOBS is expanded with file(GLOB) at configure time, so new
##    source files added to a vendored dependency are only picked up on the
##    next `cmake .` re-configure. Its objects are appended to gismo_EXTENSIONS
##    (see below) and reported back through <Name>_OBJECTS. TARGET is,
##    exactly as on the HEADER_ONLY branch, an `add_library(... INTERFACE
##    IMPORTED GLOBAL)` carrying the include dir - never an ALIAS: an ALIAS of
##    an OBJECT library, and target_link_libraries() on an OBJECT library,
##    both need CMake newer than this project's `3.1...3.10` floor
##    (CMakeLists.txt:9-13). Consumers still write
##    `target_link_libraries(<module> ... <Ns>::<Name>)` to get the include
##    dir; the compiled objects themselves arrive separately, through
##    gismo_EXTENSIONS. Unless EXPORT_SYMBOLS is given (see below), no
##    visibility property is set on `gismo_dep_<Name>`: it inherits whatever
##    CMAKE_CXX_VISIBILITY_PRESET/CMAKE_C_VISIBILITY_PRESET/
##    CMAKE_VISIBILITY_INLINES_HIDDEN are in effect at the call site.
##
## Why OBJECT + gismo_EXTENSIONS rather than a separate linked library:
## gismo_EXTENSIONS entries are sources of both gismo_static
## (cmake/gsLibrary.cmake:18-22) and shared gismo (:114-118), so the objects
## are archived into libgismo.a *and* linked into libgismo.so, and they are
## also sources of every executable in the GISMO_BUILD_LIB=OFF header-only
## mode (`add_executable(${FNAME} ${FILE} ${gismo_SOURCES} ${gismo_EXTENSIONS}
## ...)`, cmake/gismoUse.cmake:31). A library reached only through
## gismo_LINKER's target_link_libraries() would be neither archived into
## libgismo.a nor exported (`export(TARGETS gismo gismo_static ...)` names
## only those two targets, cmake/gsInstall.cmake:39-42; `install(EXPORT)` is
## commented out, :169).
##
## EXPORT_SYMBOLS: a no-value keyword. On a vendored MODE SOURCES dependency
## it sets, on `gismo_dep_<Name>`: CXX_VISIBILITY_PRESET default,
## C_VISIBILITY_PRESET default, VISIBILITY_INLINES_HIDDEN OFF - overriding
## the hidden-by-default preset cmake/gsConfig.cmake:18-23 sets for GNU,
## non-Darwin builds (CMP0063 is NEW, CMakeLists.txt:29-31, so that preset
## also governs OBJECT libraries). Some vendored code needs its symbols kept
## visible for consumers of shared libgismo.so to resolve against - e.g. the
## header-inline `gsHLBFGS<T>` calls the non-template `HLBFGS()`/
## `INIT_HLBFGS` defined in HLBFGS's own vendored .cpp files. Others must NOT
## be exported this way - vendored zlib/gzstream symbols exported from
## libgismo.so can clash with a consumer's own zlib. Combined with
## MODE HEADER_ONLY, EXPORT_SYMBOLS is a FATAL_ERROR raised in argument
## validation, before resolution or fetching: a HEADER_ONLY dependency
## compiles no code of its own, so there is no target whose visibility could
## be set, and the mismatch is visible in the call itself, independent of the
## machine it runs on - unlike whether find_package() succeeds. On a found
## dependency EXPORT_SYMBOLS is a documented no-op: a prebuilt system
## library's visibility was fixed when it was built, and find_package()
## success is machine-dependent, so erroring there would make the same
## declaration valid on one box and invalid on another.
## GISMO_INCLUDE_DIRS is appended in every case, each element deduplicated
## individually - a resolved include dir or library list can itself hold
## several paths, and the body re-runs on every `cmake .` re-configure - so
## repeated calls/re-configures cannot grow the list. gismo_EXTENSIONS
## receives $<TARGET_OBJECTS:gismo_dep_<Name>> only for a vendored MODE
## SOURCES dependency, deduplicated the same way. gismo_LINKER (the legacy
## single-library build's global, consumed via target_link_libraries() of
## `gismo`/`gismo_static`, cmake/gsLibrary.cmake:27,179) receives only a found
## (non-vendored) package's reported libraries - a vendored dependency, of
## either MODE, contributes nothing to gismo_LINKER.
##
## GIT_REPOSITORY and URL must be https:// - CI has no SSH credentials.
##
## MODE EXTERNAL_PROJECT is not supported yet: gsIpOpt, gsTrilinos and
## gsOpenCascade still carry their own ExternalProject-based fetch code.
##
## TRAP: call this only after CMakeLists.txt's `set(GISMO_INCLUDE_DIRS ...)`
## (currently line 216), which unconditionally overwrites the variable. A
## call placed between `include(gsFetch)` (currently line 191) and that
## `set()` silently loses its include-dir contribution to GISMO_INCLUDE_DIRS
## while keeping its gismo_LINKER / gismo_EXTENSIONS contribution, since the
## globals are appended independently and only GISMO_INCLUDE_DIRS is later
## clobbered.
##
## Examples:
##   gismo_add_dependency(Spectra
##     FIND_PACKAGE   spectra
##     GIT_REPOSITORY https://github.com/yixuan/spectra.git
##     GIT_TAG        v1.2.0
##     MODE           HEADER_ONLY
##     INCLUDE_SUBDIR include
##     TARGET         Spectra::Spectra)
##
##   gismo_add_dependency(HLBFGS
##     SOURCE_DIR     "${gismo_SOURCE_DIR}/external"
##     MODE           SOURCES
##     SOURCE_GLOBS   HLBFGS/HLBFGS.cpp HLBFGS/HLBFGS_BLAS.cpp HLBFGS/ICFS.cpp HLBFGS/LineSearch.cpp
##     EXPORT_SYMBOLS
##     TARGET         HLBFGS::HLBFGS)
##
function(gismo_add_dependency GAD_NAME)
  # Re-entrancy latch: a GLOBAL property (not a cache variable) so that a
  # second `cmake .` re-configure still re-runs this function - cache
  # variables persist across configures but targets do not, so latching on
  # a cache variable would leave a later configure without its target.
  get_property(_gad_resolved GLOBAL PROPERTY gismo_dependency_${GAD_NAME}_resolved SET)
  if(_gad_resolved)
    return()
  endif()

  cmake_parse_arguments(GAD "EXPORT_SYMBOLS"
    "MODE;INCLUDE_SUBDIR;TARGET;GIT_REPOSITORY;GIT_TAG;URL;URL_HASH;SOURCE_DIR"
    "FIND_PACKAGE;SOURCE_GLOBS" ${ARGN})

  # --- argument validation --------------------------------------------------
  if(NOT DEFINED GAD_MODE OR GAD_MODE STREQUAL "")
    message(FATAL_ERROR "gismo_add_dependency(${GAD_NAME}): MODE is required (HEADER_ONLY or SOURCES)")
  endif()
  if(GAD_MODE STREQUAL "EXTERNAL_PROJECT")
    message(FATAL_ERROR "gismo_add_dependency(${GAD_NAME}): MODE EXTERNAL_PROJECT is not supported yet; gsIpOpt, gsTrilinos and gsOpenCascade still carry their own ExternalProject code")
  endif()
  if(NOT GAD_MODE STREQUAL "HEADER_ONLY" AND NOT GAD_MODE STREQUAL "SOURCES")
    message(FATAL_ERROR "gismo_add_dependency(${GAD_NAME}): MODE must be HEADER_ONLY or SOURCES (got \"${GAD_MODE}\")")
  endif()
  if(GAD_EXPORT_SYMBOLS AND GAD_MODE STREQUAL "HEADER_ONLY")
    message(FATAL_ERROR "gismo_add_dependency(${GAD_NAME}): EXPORT_SYMBOLS applies only to MODE SOURCES; a MODE HEADER_ONLY dependency compiles no code whose symbol visibility could be set")
  endif()
  if(NOT DEFINED GAD_TARGET OR GAD_TARGET STREQUAL "")
    message(FATAL_ERROR "gismo_add_dependency(${GAD_NAME}): TARGET is required")
  endif()
  if(GAD_MODE STREQUAL "SOURCES" AND NOT GAD_SOURCE_GLOBS)
    message(FATAL_ERROR "gismo_add_dependency(${GAD_NAME}): MODE SOURCES requires SOURCE_GLOBS")
  endif()

  set(_gad_origin_count 0)
  if(DEFINED GAD_GIT_REPOSITORY AND NOT GAD_GIT_REPOSITORY STREQUAL "")
    math(EXPR _gad_origin_count "${_gad_origin_count}+1")
  endif()
  if(DEFINED GAD_URL AND NOT GAD_URL STREQUAL "")
    math(EXPR _gad_origin_count "${_gad_origin_count}+1")
  endif()
  if(DEFINED GAD_SOURCE_DIR AND NOT GAD_SOURCE_DIR STREQUAL "")
    math(EXPR _gad_origin_count "${_gad_origin_count}+1")
  endif()
  if(_gad_origin_count GREATER 1)
    message(FATAL_ERROR "gismo_add_dependency(${GAD_NAME}): give at most one of GIT_REPOSITORY, URL, SOURCE_DIR")
  endif()

  if(DEFINED GAD_GIT_REPOSITORY AND NOT GAD_GIT_REPOSITORY STREQUAL "" AND NOT GAD_GIT_REPOSITORY MATCHES "^https://")
    message(FATAL_ERROR "gismo_add_dependency(${GAD_NAME}): GIT_REPOSITORY must be a https:// URL - CI has no SSH credentials (got \"${GAD_GIT_REPOSITORY}\")")
  endif()
  if(DEFINED GAD_URL AND NOT GAD_URL STREQUAL "" AND NOT GAD_URL MATCHES "^https://")
    message(FATAL_ERROR "gismo_add_dependency(${GAD_NAME}): URL must be a https:// URL - CI has no SSH credentials (got \"${GAD_URL}\")")
  endif()
  if(DEFINED GAD_SOURCE_DIR AND NOT GAD_SOURCE_DIR STREQUAL "" AND NOT EXISTS "${GAD_SOURCE_DIR}")
    message(FATAL_ERROR "gismo_add_dependency(${GAD_NAME}): SOURCE_DIR \"${GAD_SOURCE_DIR}\" does not exist")
  endif()

  # --- effective fetch mode --------------------------------------------------
  string(TOUPPER "${GAD_NAME}" _gad_name_upper)
  set(_gad_override_var "GISMO_${_gad_name_upper}_FETCH")
  if(DEFINED ${_gad_override_var} AND NOT "${${_gad_override_var}}" STREQUAL "")
    set(_gad_fetch_mode "${${_gad_override_var}}")
  elseif(DEFINED GISMO_DEPENDENCY_FETCH AND NOT GISMO_DEPENDENCY_FETCH STREQUAL "")
    set(_gad_fetch_mode "${GISMO_DEPENDENCY_FETCH}")
  else()
    set(_gad_fetch_mode "AUTO")
  endif()
  if(NOT _gad_fetch_mode STREQUAL "AUTO" AND NOT _gad_fetch_mode STREQUAL "ALWAYS" AND NOT _gad_fetch_mode STREQUAL "NEVER")
    message(FATAL_ERROR "gismo_add_dependency(${GAD_NAME}): ${_gad_override_var}/GISMO_DEPENDENCY_FETCH must be one of AUTO, ALWAYS, NEVER (got \"${_gad_fetch_mode}\")")
  endif()

  set(_gad_found FALSE)
  set(_gad_vendored FALSE)
  set(_gad_include_dir "")
  set(_gad_libs "")
  set(_gad_objects "")

  # --- find_package branch ---------------------------------------------------
  if(DEFINED GAD_FIND_PACKAGE AND NOT _gad_fetch_mode STREQUAL "ALWAYS")
    list(GET GAD_FIND_PACKAGE 0 _gad_pkg)
    string(TOUPPER "${_gad_pkg}" _gad_pkg_upper)

    find_package(${GAD_FIND_PACKAGE} QUIET)

    set(_gad_found_var1 "${_gad_pkg}_FOUND")
    set(_gad_found_var2 "${_gad_pkg_upper}_FOUND")
    if(${_gad_found_var1} OR ${_gad_found_var2})

      if(TARGET ${GAD_TARGET})
        get_target_property(_gad_include_dir ${GAD_TARGET} INTERFACE_INCLUDE_DIRECTORIES)
        if(_gad_include_dir STREQUAL "_gad_include_dir-NOTFOUND")
          set(_gad_include_dir "")
        endif()
      else()
        foreach(_gad_cand
            "${_gad_pkg_upper}_INCLUDE_DIR" "${_gad_pkg}_INCLUDE_DIR"
            "${_gad_pkg_upper}_INCLUDE_DIRS" "${_gad_pkg}_INCLUDE_DIRS")
          if(_gad_include_dir STREQUAL "" AND NOT "${${_gad_cand}}" STREQUAL "")
            set(_gad_include_dir "${${_gad_cand}}")
          endif()
        endforeach()
        if(_gad_include_dir STREQUAL "")
          message(FATAL_ERROR "gismo_add_dependency(${GAD_NAME}): package ${_gad_pkg} was found but no include directory could be determined")
        endif()
      endif()

      # A found package's own CMake config routinely exports
      # INTERFACE_INCLUDE_DIRECTORIES with generator expressions the *consumer*
      # (this project) cannot evaluate at configure time - only the generator
      # can, once the build type / install prefix are known. Normalise what can
      # be normalised now instead of letting EXISTS below choke on it:
      #  - $<BUILD_INTERFACE:X>  -> X (the build-tree path a consumer here needs)
      #  - $<INSTALL_INTERFACE:X> -> dropped (relative to the OTHER project's
      #    install prefix, meaningless in this build)
      #  - any other $<...> -> kept verbatim (it may still resolve correctly at
      #    generate time) and exempted from the EXISTS check, since a generator
      #    expression is never a path that exists on disk today.
      #
      # Two standing hazards a naive per-element or whole-string match cannot
      # handle correctly, both of which real package configs produce:
      #  - a genex argument containing a semicolon, e.g. $<BUILD_INTERFACE:/a;/b>,
      #    is split by CMake's own list-splitting into the two list elements
      #    "$<BUILD_INTERFACE:/a" and "/b>" before this property value ever
      #    reaches here - neither fragment is meaningful on its own, so a split
      #    argument must be rejoined before it can be parsed.
      #  - several genexes concatenated in one element, e.g.
      #    $<BUILD_INTERFACE:/exists>$<INSTALL_INTERFACE:include> (produced by
      #    a single target_include_directories() call with both BUILD_INTERFACE
      #    and INSTALL_INTERFACE arguments), so one element can hold more than
      #    one token and each must be located and classified on its own.
      # A genex's own argument can itself be an arbitrarily nested $<...> (e.g.
      # a BUILD_INTERFACE wrapping a $<IF:...>), so locating a token's end by
      # counting $</> nesting depth is required - a non-nesting-aware stand-in
      # such as [^>]* would stop at the first '>', inside the nested expression.
      # The rejoin below tracks depth across list elements (the accumulator is
      # stored with its semicolons backslash-escaped, since list(APPEND) would
      # otherwise re-split it exactly like the original property value was);
      # the token walk further down tracks depth within one (possibly rejoined)
      # element to find each token's own matching '>'.
      set(_gad_include_dir_raw "${_gad_include_dir}")
      set(_gad_include_dir "")
      set(_gad_genex_resolved FALSE)

      # --- rejoin list elements a semicolon-bearing genex argument was split into ---
      set(_gad_joined_elements "")
      set(_gad_pending "")
      foreach(_gad_raw_dir ${_gad_include_dir_raw})
        if(_gad_pending STREQUAL "")
          set(_gad_candidate "${_gad_raw_dir}")
        else()
          set(_gad_candidate "${_gad_pending};${_gad_raw_dir}")
        endif()
        string(REGEX MATCHALL "\\$<" _gad_opens "${_gad_candidate}")
        list(LENGTH _gad_opens _gad_nopen)
        string(REGEX MATCHALL ">" _gad_closes "${_gad_candidate}")
        list(LENGTH _gad_closes _gad_nclose)
        if(_gad_nopen GREATER _gad_nclose)
          # more opens than closes seen so far - the genex is still split
          # across a later list element; keep accumulating.
          set(_gad_pending "${_gad_candidate}")
        else()
          string(REPLACE ";" "\;" _gad_candidate_esc "${_gad_candidate}")
          list(APPEND _gad_joined_elements "${_gad_candidate_esc}")
          set(_gad_pending "")
        endif()
      endforeach()
      if(NOT _gad_pending STREQUAL "")
        # unbalanced to the end (malformed genex) - keep it as-is; the
        # per-element parse below will fail to find a matching '>' and falls
        # through to "kept, exempt, reported" for this element too.
        string(REPLACE ";" "\;" _gad_pending_esc "${_gad_pending}")
        list(APPEND _gad_joined_elements "${_gad_pending_esc}")
      endif()

      # --- walk the genex tokens of each (possibly rejoined) element ---
      foreach(_gad_elem ${_gad_joined_elements})
        set(_gad_norm_elem "")
        set(_gad_remaining "${_gad_elem}")
        while(NOT _gad_remaining STREQUAL "")
          string(FIND "${_gad_remaining}" "$<" _gad_open_pos)
          if(_gad_open_pos EQUAL -1)
            set(_gad_norm_elem "${_gad_norm_elem}${_gad_remaining}")
            set(_gad_remaining "")
          else()
            if(_gad_open_pos GREATER 0)
              string(SUBSTRING "${_gad_remaining}" 0 ${_gad_open_pos} _gad_lit)
              set(_gad_norm_elem "${_gad_norm_elem}${_gad_lit}")
            endif()
            # find this token's own matching '>' by nesting depth, not by regex
            string(LENGTH "${_gad_remaining}" _gad_len)
            set(_gad_i ${_gad_open_pos})
            set(_gad_depth 0)
            set(_gad_tok_end -1)
            while(_gad_i LESS _gad_len)
              string(SUBSTRING "${_gad_remaining}" ${_gad_i} 2 _gad_two)
              if(_gad_two STREQUAL "$<")
                math(EXPR _gad_depth "${_gad_depth}+1")
                math(EXPR _gad_i "${_gad_i}+2")
              else()
                string(SUBSTRING "${_gad_remaining}" ${_gad_i} 1 _gad_ch)
                math(EXPR _gad_i "${_gad_i}+1")
                if(_gad_ch STREQUAL ">")
                  math(EXPR _gad_depth "${_gad_depth}-1")
                  if(_gad_depth EQUAL 0)
                    set(_gad_tok_end ${_gad_i})
                    break()
                  endif()
                endif()
              endif()
            endwhile()
            if(_gad_tok_end EQUAL -1)
              # no matching '>' - malformed genex; still not a path EXISTS can
              # check, so keep and report it exactly like a bare unrecognised
              # genex, rather than letting it through unannounced.
              string(SUBSTRING "${_gad_remaining}" ${_gad_open_pos} -1 _gad_tok)
              message(STATUS "gismo_add_dependency(${GAD_NAME}): include directory \"${_gad_tok}\" is an unevaluated generator expression - kept on ${GAD_TARGET}, exempt from the EXISTS check below")
              set(_gad_norm_elem "${_gad_norm_elem}${_gad_tok}")
              set(_gad_remaining "")
            else()
              math(EXPR _gad_tok_len "${_gad_tok_end}-${_gad_open_pos}")
              string(SUBSTRING "${_gad_remaining}" ${_gad_open_pos} ${_gad_tok_len} _gad_tok)
              string(SUBSTRING "${_gad_remaining}" ${_gad_tok_end} -1 _gad_remaining)

              if(_gad_tok MATCHES "^\\$<BUILD_INTERFACE:(.*)>$")
                set(_gad_dir_payload "${CMAKE_MATCH_1}")
                if(_gad_dir_payload MATCHES "\\$<")
                  # the BUILD_INTERFACE payload is itself an unresolved generator
                  # expression (e.g. a nested $<IF:...>) - still not a path EXISTS
                  # can check, so report it exactly like the bare-genex case below.
                  message(STATUS "gismo_add_dependency(${GAD_NAME}): include directory \"${_gad_dir_payload}\" is an unevaluated generator expression - kept on ${GAD_TARGET}, exempt from the EXISTS check below")
                endif()
                set(_gad_norm_elem "${_gad_norm_elem}${_gad_dir_payload}")
                set(_gad_genex_resolved TRUE)
              elseif(_gad_tok MATCHES "^\\$<INSTALL_INTERFACE:.*>$")
                set(_gad_genex_resolved TRUE)
              else()
                message(STATUS "gismo_add_dependency(${GAD_NAME}): include directory \"${_gad_tok}\" is an unevaluated generator expression - kept on ${GAD_TARGET}, exempt from the EXISTS check below")
                set(_gad_norm_elem "${_gad_norm_elem}${_gad_tok}")
              endif()
            endif()
          endif()
        endwhile()
        # a normalised element that collapsed to nothing (a pure, dropped
        # $<INSTALL_INTERFACE:...>) must not become an empty list entry -
        # list(APPEND var "") appends a real (empty-string) element.
        if(NOT _gad_norm_elem STREQUAL "")
          list(APPEND _gad_include_dir "${_gad_norm_elem}")
        endif()
      endforeach()

      if(_gad_genex_resolved AND TARGET ${GAD_TARGET})
        # the target already existed (defined by the package's own config) -
        # push the normalised list back onto it, it is not rewritten below.
        set_target_properties(${GAD_TARGET} PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${_gad_include_dir}")
      endif()

      foreach(_gad_cand
          "${_gad_pkg_upper}_LIBRARY" "${_gad_pkg}_LIBRARY"
          "${_gad_pkg_upper}_LIBRARIES" "${_gad_pkg}_LIBRARIES")
        if(_gad_libs STREQUAL "" AND NOT "${${_gad_cand}}" STREQUAL "")
          set(_gad_libs "${${_gad_cand}}")
        endif()
      endforeach()

      # A found package whose Config.cmake defines its imported target
      # directly (add_library(...IMPORTED) + IMPORTED_LOCATION, the shape
      # CMake itself recommends) routinely reports no <Pkg>_LIBRARY/_LIBRARIES
      # variable at all - the loop above finds nothing even though TARGET
      # carries a perfectly usable link interface. Fall back to reading it
      # off TARGET itself, but only when the variable-based lookup came up
      # empty: a package that does report *_LIBRARY/_LIBRARIES is trusted on
      # that documented contract, as this function always has.
      if(_gad_libs STREQUAL "" AND TARGET ${GAD_TARGET})
        set(_gad_target_visited "")
        _gismo_dependency_link_paths("${GAD_NAME}" "${GAD_TARGET}" "" _gad_target_libs _gad_target_visited)
        if(NOT _gad_target_libs STREQUAL "")
          set(_gad_libs "${_gad_target_libs}")
        endif()
      endif()

      if(NOT TARGET ${GAD_TARGET})
        add_library(${GAD_TARGET} INTERFACE IMPORTED GLOBAL)
        set_target_properties(${GAD_TARGET} PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${_gad_include_dir}")
        if(NOT _gad_libs STREQUAL "")
          set_target_properties(${GAD_TARGET} PROPERTIES INTERFACE_LINK_LIBRARIES "${_gad_libs}")
        endif()
      endif()

      if(_gad_include_dir STREQUAL "")
        message(FATAL_ERROR "gismo_add_dependency(${GAD_NAME}): package ${_gad_pkg} was found but no include directory could be determined")
      endif()
      foreach(_gad_dir ${_gad_include_dir})
        # an unresolved generator expression is never a path that exists on
        # disk at configure time (reported above instead, when it first
        # survived normalisation) - only plain paths are checked here.
        if(NOT _gad_dir MATCHES "\\$<" AND NOT EXISTS "${_gad_dir}")
          message(FATAL_ERROR "gismo_add_dependency(${GAD_NAME}): include directory \"${_gad_dir}\" does not exist")
        endif()
      endforeach()

      set(_gad_vendored FALSE)
      set(_gad_found TRUE)
      message(STATUS "Using ${GAD_NAME} at ${_gad_include_dir}")
    endif()
  endif()

  # --- fetch/vendored branch --------------------------------------------------
  if(NOT _gad_found)
    if(DEFINED GAD_SOURCE_DIR AND NOT GAD_SOURCE_DIR STREQUAL "")
      set(_gad_source_dir "${GAD_SOURCE_DIR}")
    else()
      if(_gad_fetch_mode STREQUAL "NEVER")
        message(FATAL_ERROR "gismo_add_dependency(${GAD_NAME}): dependency was not found and fetching is disabled (${_gad_override_var}/GISMO_DEPENDENCY_FETCH = NEVER); point <pkg>_DIR/CMAKE_PREFIX_PATH at an installed copy, or allow fetching")
      endif()

      if(DEFINED GAD_GIT_REPOSITORY AND NOT GAD_GIT_REPOSITORY STREQUAL "")
        set(_gad_origin_args GIT_REPOSITORY "${GAD_GIT_REPOSITORY}")
        if(DEFINED GAD_GIT_TAG AND NOT GAD_GIT_TAG STREQUAL "")
          list(APPEND _gad_origin_args GIT_TAG "${GAD_GIT_TAG}")
        endif()
      elseif(DEFINED GAD_URL AND NOT GAD_URL STREQUAL "")
        set(_gad_origin_args URL "${GAD_URL}")
        if(DEFINED GAD_URL_HASH AND NOT GAD_URL_HASH STREQUAL "")
          list(APPEND _gad_origin_args URL_HASH "${GAD_URL_HASH}")
        endif()
      else()
        message(FATAL_ERROR "gismo_add_dependency(${GAD_NAME}): no GIT_REPOSITORY, URL or SOURCE_DIR given and no FIND_PACKAGE succeeded")
      endif()

      gismo_fetch_directory(${GAD_NAME} ${_gad_origin_args} DESTINATION external)
      set(_gad_source_dir "${gismo_SOURCE_DIR}/external/${GAD_NAME}")
    endif()

    if(DEFINED GAD_INCLUDE_SUBDIR AND NOT GAD_INCLUDE_SUBDIR STREQUAL "")
      set(_gad_include_dir "${_gad_source_dir}/${GAD_INCLUDE_SUBDIR}")
    else()
      set(_gad_include_dir "${_gad_source_dir}")
    endif()

    if(NOT EXISTS "${_gad_include_dir}")
      message(FATAL_ERROR "gismo_add_dependency(${GAD_NAME}): include directory \"${_gad_include_dir}\" does not exist")
    endif()

    if(GAD_MODE STREQUAL "HEADER_ONLY")
      if(NOT TARGET ${GAD_TARGET})
        add_library(${GAD_TARGET} INTERFACE IMPORTED GLOBAL)
        set_target_properties(${GAD_TARGET} PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${_gad_include_dir}")
      endif()
    else() # SOURCES
      set(_gad_sources "")
      foreach(_gad_glob ${GAD_SOURCE_GLOBS})
        file(GLOB _gad_glob_result "${_gad_source_dir}/${_gad_glob}")
        list(APPEND _gad_sources ${_gad_glob_result})
      endforeach()
      if(NOT _gad_sources)
        message(FATAL_ERROR "gismo_add_dependency(${GAD_NAME}): SOURCE_GLOBS (${GAD_SOURCE_GLOBS}) matched no files under \"${_gad_source_dir}\"")
      endif()
      if(NOT TARGET gismo_dep_${GAD_NAME})
        add_library(gismo_dep_${GAD_NAME} OBJECT ${_gad_sources})
        target_include_directories(gismo_dep_${GAD_NAME} PRIVATE "${_gad_include_dir}")
        set_target_properties(gismo_dep_${GAD_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)
        if(MSVC)
          target_compile_options(gismo_dep_${GAD_NAME} PRIVATE /W0)
        elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
          target_compile_options(gismo_dep_${GAD_NAME} PRIVATE -w)
        endif()
        if(GAD_EXPORT_SYMBOLS)
          set_target_properties(gismo_dep_${GAD_NAME} PROPERTIES
            CXX_VISIBILITY_PRESET default
            C_VISIBILITY_PRESET default
            VISIBILITY_INLINES_HIDDEN OFF)
        endif()
      endif()
      if(NOT TARGET ${GAD_TARGET})
        add_library(${GAD_TARGET} INTERFACE IMPORTED GLOBAL)
        set_target_properties(${GAD_TARGET} PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${_gad_include_dir}")
      endif()
      set(_gad_objects "$<TARGET_OBJECTS:gismo_dep_${GAD_NAME}>")
    endif()

    set(_gad_vendored TRUE)
    message(STATUS "Using fetched ${GAD_NAME} at ${_gad_include_dir}")
  endif()

  # --- outputs -----------------------------------------------------------
  # Reaching this point means resolution succeeded: every failure path above
  # is a FATAL_ERROR, which stops the configure outright.
  set(${GAD_NAME}_FOUND TRUE CACHE INTERNAL "Whether ${GAD_NAME} was resolved by gismo_add_dependency()")
  set(${GAD_NAME}_VENDORED ${_gad_vendored} CACHE INTERNAL "Whether ${GAD_NAME} was fetched/vendored by G+Smo")
  set(${GAD_NAME}_OBJECTS "${_gad_objects}" CACHE INTERNAL "Objects of ${GAD_NAME} compiled by G+Smo (empty unless vendored MODE SOURCES)")

  # --- legacy single-library build globals --------------------------------
  # _gad_include_dir / _gad_linker_entry / _gad_objects can each be a
  # semicolon-separated list in their own right (a found package's
  # *_INCLUDE_DIRS / *_LIBRARIES, or a target's INTERFACE_INCLUDE_DIRECTORIES,
  # routinely hold more than one path) - list(FIND ...) against the whole
  # value would compare it as a single element and never match, so dedup
  # element-by-element instead. CACHE INTERNAL always overwrites regardless
  # of FORCE, so it is safe (and necessary, since the function body re-runs
  # on every `cmake .`) to just recompute and re-set the deduplicated lists
  # (GISMO_INCLUDE_DIRS, gismo_LINKER, gismo_EXTENSIONS) unconditionally.
  set(_gad_new_include_dirs "${GISMO_INCLUDE_DIRS}")
  foreach(_gad_one_dir ${_gad_include_dir})
    list(FIND _gad_new_include_dirs "${_gad_one_dir}" _gad_incdir_idx)
    if(_gad_incdir_idx EQUAL -1)
      list(APPEND _gad_new_include_dirs "${_gad_one_dir}")
    endif()
  endforeach()
  set(GISMO_INCLUDE_DIRS ${_gad_new_include_dirs} CACHE INTERNAL "Gismo include directories")

  set(_gad_linker_entry "")
  if(NOT _gad_vendored AND NOT "${_gad_libs}" STREQUAL "")
    set(_gad_linker_entry "${_gad_libs}")
  endif()
  if(NOT _gad_linker_entry STREQUAL "")
    set(_gad_new_linker "${gismo_LINKER}")
    foreach(_gad_one_lib ${_gad_linker_entry})
      list(FIND _gad_new_linker "${_gad_one_lib}" _gad_linker_idx)
      if(_gad_linker_idx EQUAL -1)
        list(APPEND _gad_new_linker "${_gad_one_lib}")
      endif()
    endforeach()
    set(gismo_LINKER ${_gad_new_linker} CACHE INTERNAL "${PROJECT_NAME} extra linker objects")
  endif()

  if(NOT _gad_objects STREQUAL "")
    set(_gad_new_extensions "${gismo_EXTENSIONS}")
    foreach(_gad_one_object ${_gad_objects})
      list(FIND _gad_new_extensions "${_gad_one_object}" _gad_extensions_idx)
      if(_gad_extensions_idx EQUAL -1)
        list(APPEND _gad_new_extensions "${_gad_one_object}")
      endif()
    endforeach()
    set(gismo_EXTENSIONS ${_gad_new_extensions} CACHE INTERNAL "Gismo extensions to be included")
  endif()

  set_property(GLOBAL PROPERTY gismo_dependency_${GAD_NAME}_resolved TRUE)
endfunction()

## _gismo_dependency_link_paths(<label> <target> <visited-in> <out-libs-var> <out-visited-var>)
##
## Resolves the absolute link-library paths of a found package's own imported
## target - the modern Config.cmake shape (add_library(...IMPORTED) +
## IMPORTED_LOCATION) that reports nothing through the legacy
## <Pkg>_LIBRARY/_LIBRARIES variables gismo_add_dependency() otherwise reads.
## Returns, in order: TARGET's own compiled-artifact location (if it has one -
## see the IMPORTED_LOCATION selection below), then every entry of its
## INTERFACE_LINK_LIBRARIES that is a usable linker item - an absolute path
## or a plain token (e.g. "pthread") kept verbatim (the JIT config's -l
## rewrite, cmake/gsJITConfigXml.cmake:112-117, is correct for a bare name as-
## is), or the name of another target (imported or not), recursed into. An
## entry containing an unevaluated generator expression cannot be resolved at
## configure time and is skipped, with one status message per skipped entry.
##
## IMPORTED_LOCATION selection, in order, each step tried only when the
## previous one produced nothing: IMPORTED_LOCATION_<CMAKE_BUILD_TYPE
## uppercased> when CMAKE_BUILD_TYPE is set and that property is populated;
## else plain IMPORTED_LOCATION when set; else the first populated
## IMPORTED_LOCATION_<CONFIG> over the target's own IMPORTED_CONFIGURATIONS.
## On a multi-config generator CMAKE_BUILD_TYPE is normally unset, which rules
## out only the first step - the second (plain IMPORTED_LOCATION) still wins
## whenever the target's Config.cmake sets it, same as on a single-config
## build; only a target with no plain IMPORTED_LOCATION at all falls through
## to the third, genuinely config-independent step. Whichever step fires, the
## result cannot vary per configuration within one configure run, since this
## function runs once at configure time, before a build type is selected. An
## INTERFACE_LIBRARY target has no compiled artifact of its own; querying
## IMPORTED_LOCATION* on one is a hard error on older CMake, so TYPE is
## checked first rather than relying on a NOTFOUND result.
##
## <visited-in> is the semicolon-separated set of every target name visited
## so far on this call's whole recursion tree (not just the current path back
## to the root) - the cycle guard for a Config.cmake whose
## INTERFACE_LINK_LIBRARIES names a target that (directly or transitively)
## links back to a target already walked, e.g. two targets A and B where A's
## INTERFACE_LINK_LIBRARIES names B and B's names A. Threading the full
## visited set (rather than only the current path) also means a diamond - two
## branches both naming the same third target - is walked once, not once per
## branch. CMake functions have no reference parameters, so like the resolved
## paths, the updated visited set is returned through an out-variable name,
## written with PARENT_SCOPE.
function(_gismo_dependency_link_paths GLP_LABEL GLP_TARGET GLP_VISITED_IN GLP_OUT_LIBS GLP_OUT_VISITED)
  set(_glp_visited "${GLP_VISITED_IN}")
  list(FIND _glp_visited "${GLP_TARGET}" _glp_visited_idx)
  if(NOT _glp_visited_idx EQUAL -1)
    set(${GLP_OUT_LIBS} "" PARENT_SCOPE)
    set(${GLP_OUT_VISITED} "${_glp_visited}" PARENT_SCOPE)
    return()
  endif()
  list(APPEND _glp_visited "${GLP_TARGET}")

  set(_glp_result "")

  get_target_property(_glp_type ${GLP_TARGET} TYPE)
  if(NOT _glp_type STREQUAL "INTERFACE_LIBRARY")
    set(_glp_location "")
    if(CMAKE_BUILD_TYPE)
      string(TOUPPER "${CMAKE_BUILD_TYPE}" _glp_bt_upper)
      get_target_property(_glp_location ${GLP_TARGET} IMPORTED_LOCATION_${_glp_bt_upper})
      if(_glp_location STREQUAL "_glp_location-NOTFOUND")
        set(_glp_location "")
      endif()
    endif()
    if(_glp_location STREQUAL "")
      get_target_property(_glp_location ${GLP_TARGET} IMPORTED_LOCATION)
      if(_glp_location STREQUAL "_glp_location-NOTFOUND")
        set(_glp_location "")
      endif()
    endif()
    if(_glp_location STREQUAL "")
      get_target_property(_glp_configs ${GLP_TARGET} IMPORTED_CONFIGURATIONS)
      if(_glp_configs STREQUAL "_glp_configs-NOTFOUND")
        set(_glp_configs "")
      endif()
      foreach(_glp_cfg ${_glp_configs})
        if(_glp_location STREQUAL "")
          get_target_property(_glp_cand ${GLP_TARGET} IMPORTED_LOCATION_${_glp_cfg})
          if(NOT _glp_cand STREQUAL "_glp_cand-NOTFOUND" AND NOT _glp_cand STREQUAL "")
            set(_glp_location "${_glp_cand}")
          endif()
        endif()
      endforeach()
    endif()
    if(NOT _glp_location STREQUAL "")
      list(APPEND _glp_result "${_glp_location}")
    endif()
  endif()

  get_target_property(_glp_iface ${GLP_TARGET} INTERFACE_LINK_LIBRARIES)
  if(_glp_iface STREQUAL "_glp_iface-NOTFOUND")
    set(_glp_iface "")
  endif()

  # A generator expression's own argument can contain a semicolon (e.g.
  # $<$<PLATFORM_ID:Linux>:/a;/b>) - CMake's list-splitting breaks this into
  # separate list elements ("$<$<PLATFORM_ID:Linux>:/a" and "/b>") before this
  # property value ever reaches here, and neither fragment is a meaningful
  # linker item on its own: the first still looks like an open genex (kept,
  # correctly), but the second no longer contains "$<" and would otherwise be
  # misclassified as a plain absolute path/token and leak into gismo_LINKER.
  # Rejoin first - same idiom already used for include dirs above (the
  # "rejoin list elements a semicolon-bearing genex argument was split into"
  # block): accumulate consecutive elements while "$<" count exceeds ">"
  # count, then store the joined candidate with its semicolons backslash-
  # escaped, since list(APPEND) would otherwise re-split it exactly like the
  # original property value was.
  set(_glp_joined_iface "")
  set(_glp_pending "")
  foreach(_glp_raw_elem ${_glp_iface})
    if(_glp_pending STREQUAL "")
      set(_glp_candidate "${_glp_raw_elem}")
    else()
      set(_glp_candidate "${_glp_pending};${_glp_raw_elem}")
    endif()
    string(REGEX MATCHALL "\\$<" _glp_opens "${_glp_candidate}")
    list(LENGTH _glp_opens _glp_nopen)
    string(REGEX MATCHALL ">" _glp_closes "${_glp_candidate}")
    list(LENGTH _glp_closes _glp_nclose)
    if(_glp_nopen GREATER _glp_nclose)
      # more opens than closes seen so far - the genex is still split across
      # a later list element; keep accumulating.
      set(_glp_pending "${_glp_candidate}")
    else()
      string(REPLACE ";" "\;" _glp_candidate_esc "${_glp_candidate}")
      list(APPEND _glp_joined_iface "${_glp_candidate_esc}")
      set(_glp_pending "")
    endif()
  endforeach()
  if(NOT _glp_pending STREQUAL "")
    # unbalanced to the end (malformed genex) - still not a linker item this
    # function can resolve, so keep it together and let the classification
    # below skip (and report) it as one entry, rather than leaking a bare
    # trailing fragment.
    string(REPLACE ";" "\;" _glp_pending_esc "${_glp_pending}")
    list(APPEND _glp_joined_iface "${_glp_pending_esc}")
  endif()

  foreach(_glp_elem ${_glp_joined_iface})
    if(_glp_elem MATCHES "\\$<")
      message(STATUS "gismo_add_dependency(${GLP_LABEL}): INTERFACE_LINK_LIBRARIES entry \"${_glp_elem}\" on ${GLP_TARGET} is an unevaluated generator expression - skipped, not appended to gismo_LINKER")
    elseif(IS_ABSOLUTE "${_glp_elem}")
      list(APPEND _glp_result "${_glp_elem}")
    elseif(TARGET ${_glp_elem})
      _gismo_dependency_link_paths("${GLP_LABEL}" "${_glp_elem}" "${_glp_visited}" _glp_nested_libs _glp_visited)
      if(NOT _glp_nested_libs STREQUAL "")
        list(APPEND _glp_result ${_glp_nested_libs})
      endif()
    else()
      list(APPEND _glp_result "${_glp_elem}")
    endif()
  endforeach()

  set(${GLP_OUT_LIBS} "${_glp_result}" PARENT_SCOPE)
  set(${GLP_OUT_VISITED} "${_glp_visited}" PARENT_SCOPE)
endfunction()

#function(gismo_add_plugin PLUGIN)
