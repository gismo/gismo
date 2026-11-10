# - Try to find Eigen (>= 5.0), either via a CMake package config or by
#   locating its headers directly.
#
#   find_package(Eigen 5.0)
#
# Once done this will define
#
#  EIGEN_FOUND        - system has Eigen with the requested version
#  EIGEN_INCLUDE_DIR   - the Eigen include directory (the one containing Eigen/Core)
#  EIGEN_VERSION       - the Eigen version that was found
#
# Hints:
#  Eigen_DIR - path to an Eigen root (the directory containing Eigen/Core),
#              checked before any system search.
#
# Note: Eigen 5.0 dropped the old 3-level EIGEN_WORLD_VERSION/EIGEN_MAJOR_VERSION/
# EIGEN_MINOR_VERSION scheme (used by the legacy FindEigen3.cmake in this same
# directory) in favor of EIGEN_MAJOR_VERSION/EIGEN_MINOR_VERSION/EIGEN_PATCH_VERSION.

if(NOT Eigen_FIND_VERSION)
  if(NOT Eigen_FIND_VERSION_MAJOR)
    set(Eigen_FIND_VERSION_MAJOR 5)
  endif()
  if(NOT Eigen_FIND_VERSION_MINOR)
    set(Eigen_FIND_VERSION_MINOR 0)
  endif()
  if(NOT Eigen_FIND_VERSION_PATCH)
    set(Eigen_FIND_VERSION_PATCH 0)
  endif()
  set(Eigen_FIND_VERSION "${Eigen_FIND_VERSION_MAJOR}.${Eigen_FIND_VERSION_MINOR}.${Eigen_FIND_VERSION_PATCH}")
endif()

# 1) A system package with a modern CMake config (Eigen3Config.cmake, still the
#    package name used upstream even at version 5.x). Only tried when Eigen_DIR
#    is NOT set: Eigen_DIR is documented as a plain Eigen root (a directory
#    containing Eigen/Core, e.g. a raw source checkout), not a CMake install
#    prefix, and pointing find_package's CONFIG mode at one can pick up a
#    leftover/incomplete Eigen3Config.cmake (e.g. one written by configuring
#    Eigen's own build in-tree) whose include() of Eigen3Targets.cmake hard-
#    errors instead of failing gracefully.
if(NOT EIGEN_INCLUDE_DIR AND NOT Eigen_DIR)
  find_package(Eigen3 ${Eigen_FIND_VERSION} CONFIG QUIET)
  if(Eigen3_FOUND AND TARGET Eigen3::Eigen)
    get_target_property(_eigen_config_dirs Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
    # INTERFACE_INCLUDE_DIRECTORIES is a list and, for a build-tree export,
    # typically carries $<BUILD_INTERFACE:...>/$<INSTALL_INTERFACE:...>
    # generator expressions rather than plain paths. Those are evaluated at
    # CMake's generate step, which runs after this find-module, so they must
    # be resolved by hand here: unwrap BUILD_INTERFACE, skip INSTALL_INTERFACE
    # (its path is relative to an install prefix that may not exist yet) and
    # any entry still carrying "$<" after that, and take the first surviving
    # entry that is an actual Eigen root (mirrors how _eigen_check_version()
    # below recognizes one). EIGEN_INCLUDE_DIR/EIGEN_VERSION_OK are left
    # untouched when nothing qualifies, so the manual header search further
    # down runs as the fallback.
    if(_eigen_config_dirs)
      foreach(_eigen_config_dir IN LISTS _eigen_config_dirs)
        if(_eigen_config_dir MATCHES "^\\$<BUILD_INTERFACE:(.*)>$")
          set(_eigen_config_dir "${CMAKE_MATCH_1}")
        endif()
        if(_eigen_config_dir MATCHES "\\$<")
          continue()
        endif()
        if(EXISTS "${_eigen_config_dir}/Eigen/Core")
          set(EIGEN_INCLUDE_DIR "${_eigen_config_dir}")
          set(EIGEN_VERSION ${Eigen3_VERSION})
          set(EIGEN_VERSION_OK TRUE)
          break()
        endif()
      endforeach()
    endif()
    unset(_eigen_config_dirs)
    unset(_eigen_config_dir)
  endif()
endif()

# 2) Manual header search (Eigen_DIR hint, then system paths), mirroring the
#    legacy FindEigen3.cmake signature-file heuristic.
macro(_eigen_check_version)
  # The two version schemes use the same macro names for different fields:
  #   Eigen >= 5.0: Eigen/Version defines MAJOR.MINOR.PATCH (5.0.0), plus a
  #                 legacy EIGEN_WORLD_VERSION 3 that must be ignored;
  #   Eigen 3.x:    Eigen/src/Core/util/Macros.h defines WORLD.MAJOR.MINOR (3.4.0)
  #                 and there is no Eigen/Version file.
  # A pre-5.0 Eigen (e.g. a distribution's eigen3 package) must be rejected by
  # the version check below, not abort the configure by reading a missing file.
  set(EIGEN_VERSION "")
  if(EXISTS "${EIGEN_INCLUDE_DIR}/Eigen/Version")
    file(READ "${EIGEN_INCLUDE_DIR}/Eigen/Version" _eigen_version_header)
    string(REGEX MATCH "define[ \t]+EIGEN_MAJOR_VERSION[ \t]+([0-9]+)" _eigen_match "${_eigen_version_header}")
    set(EIGEN_MAJOR_VERSION "${CMAKE_MATCH_1}")
    string(REGEX MATCH "define[ \t]+EIGEN_MINOR_VERSION[ \t]+([0-9]+)" _eigen_match "${_eigen_version_header}")
    set(EIGEN_MINOR_VERSION "${CMAKE_MATCH_1}")
    string(REGEX MATCH "define[ \t]+EIGEN_PATCH_VERSION[ \t]+([0-9]+)" _eigen_match "${_eigen_version_header}")
    set(EIGEN_PATCH_VERSION "${CMAKE_MATCH_1}")
    set(EIGEN_VERSION "${EIGEN_MAJOR_VERSION}.${EIGEN_MINOR_VERSION}.${EIGEN_PATCH_VERSION}")
  elseif(EXISTS "${EIGEN_INCLUDE_DIR}/Eigen/src/Core/util/Macros.h")
    file(READ "${EIGEN_INCLUDE_DIR}/Eigen/src/Core/util/Macros.h" _eigen_version_header)
    string(REGEX MATCH "define[ \t]+EIGEN_WORLD_VERSION[ \t]+([0-9]+)" _eigen_match "${_eigen_version_header}")
    set(_eigen_world "${CMAKE_MATCH_1}")
    string(REGEX MATCH "define[ \t]+EIGEN_MAJOR_VERSION[ \t]+([0-9]+)" _eigen_match "${_eigen_version_header}")
    set(_eigen_major "${CMAKE_MATCH_1}")
    string(REGEX MATCH "define[ \t]+EIGEN_MINOR_VERSION[ \t]+([0-9]+)" _eigen_match "${_eigen_version_header}")
    set(_eigen_minor "${CMAKE_MATCH_1}")
    set(EIGEN_VERSION "${_eigen_world}.${_eigen_major}.${_eigen_minor}")
  endif()

  if(NOT EIGEN_VERSION MATCHES "^[0-9]+\\.[0-9]+\\.[0-9]+$")
    set(EIGEN_VERSION_OK FALSE)
    message(STATUS "Could not determine the Eigen version in ${EIGEN_INCLUDE_DIR}; "
                    "at least version ${Eigen_FIND_VERSION} is required")
  elseif(EIGEN_VERSION VERSION_LESS Eigen_FIND_VERSION)
    set(EIGEN_VERSION_OK FALSE)
    message(STATUS "Eigen version ${EIGEN_VERSION} found in ${EIGEN_INCLUDE_DIR}, "
                    "but at least version ${Eigen_FIND_VERSION} is required")
  else()
    set(EIGEN_VERSION_OK TRUE)
  endif()
endmacro()

if(NOT EIGEN_INCLUDE_DIR)
  find_path(EIGEN_INCLUDE_DIR NAMES signature_of_eigen3_matrix_library
    HINTS "${Eigen_DIR}"
    PATHS
    ${CMAKE_INSTALL_PREFIX}/include
    PATH_SUFFIXES eigen3 eigen
  )
endif()

# EIGEN_INCLUDE_DIR may already be set here without EIGEN_VERSION_OK having
# been checked: it is a find_path CACHE variable (survives a re-configure and
# a user-supplied -DEIGEN_INCLUDE_DIR=...), while EIGEN_VERSION_OK is a plain
# variable that a fresh CMake process does not see. NOT EIGEN_VERSION_OK also
# guards the config-package hit above (1), which already set it TRUE together
# with a single resolved EIGEN_INCLUDE_DIR, so re-deriving the version from
# headers here would just be redundant.
if(EIGEN_INCLUDE_DIR AND NOT EIGEN_VERSION_OK)
  _eigen_check_version()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Eigen DEFAULT_MSG EIGEN_INCLUDE_DIR EIGEN_VERSION_OK)

mark_as_advanced(EIGEN_INCLUDE_DIR)
