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
    get_target_property(EIGEN_INCLUDE_DIR Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
    set(EIGEN_VERSION ${Eigen3_VERSION})
    set(EIGEN_VERSION_OK TRUE)
  endif()
endif()

# 2) Manual header search (Eigen_DIR hint, then system paths), mirroring the
#    legacy FindEigen3.cmake signature-file heuristic.
macro(_eigen_check_version)
  # As of Eigen 5.0, the version macros live in Eigen/Version (included by
  # Eigen/Core), not in Eigen/src/Core/util/Macros.h.
  file(READ "${EIGEN_INCLUDE_DIR}/Eigen/Version" _eigen_version_header)

  string(REGEX MATCH "define[ \t]+EIGEN_MAJOR_VERSION[ \t]+([0-9]+)" _eigen_major_version_match "${_eigen_version_header}")
  set(EIGEN_MAJOR_VERSION "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+EIGEN_MINOR_VERSION[ \t]+([0-9]+)" _eigen_minor_version_match "${_eigen_version_header}")
  set(EIGEN_MINOR_VERSION "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+EIGEN_PATCH_VERSION[ \t]+([0-9]+)" _eigen_patch_version_match "${_eigen_version_header}")
  set(EIGEN_PATCH_VERSION "${CMAKE_MATCH_1}")

  set(EIGEN_VERSION ${EIGEN_MAJOR_VERSION}.${EIGEN_MINOR_VERSION}.${EIGEN_PATCH_VERSION})
  if(${EIGEN_VERSION} VERSION_LESS ${Eigen_FIND_VERSION})
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
  if(EIGEN_INCLUDE_DIR)
    _eigen_check_version()
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Eigen DEFAULT_MSG EIGEN_INCLUDE_DIR EIGEN_VERSION_OK)

mark_as_advanced(EIGEN_INCLUDE_DIR)
