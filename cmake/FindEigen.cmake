######################################################################
## FindEigen.cmake
## This file is part of the G+Smo library.
##
## Minimal find module for an external Eigen 5.x installation or source tree.
##
## Input:
##   Eigen_DIR    - Eigen source root, install prefix, or include directory.
##                  When set, it is the ONLY place searched: CMAKE_PREFIX_PATH,
##                  CMAKE_INCLUDE_PATH and the system paths are skipped, so an
##                  explicit choice cannot be overridden by another Eigen found
##                  elsewhere. When unset, the normal find_path search (system
##                  paths, CMAKE_PREFIX_PATH, etc.) runs.
## Output:
##   EIGEN_INCLUDE_DIR - cache path to the directory containing Eigen/ and
##                        unsupported/ (or their install-tree equivalent).
##   Eigen_VERSION     - "MAJOR.MINOR.PATCH", re-read on every configure.
##   Eigen_FOUND       - set by find_package_handle_standard_args.
##
## Only Eigen 5's own Eigen/Version header is parsed for the version: Eigen
## versions before 5 ship no such header, so a 3.4 tree is correctly reported
## as version-less rather than misread from Macros.h (whose WORLD/MAJOR
## numbering does not match the MAJOR.MINOR.PATCH scheme used from 5.0 on).
######################################################################

include(FindPackageHandleStandardArgs)
unset(Eigen_VERSION)                      # plain variable, NEVER cached: re-read on every configure
set(Eigen_REJECTED_DIR "")                # plain variable, NEVER stale: last rejected EIGEN_INCLUDE_DIR, for callers' error messages

# find_path never re-searches while its cache entry is already set, so a
# changed Eigen_DIR hint would otherwise be silently ignored on reconfigure.
if(DEFINED _Eigen_DIR_SEARCHED AND NOT "${Eigen_DIR}" STREQUAL "${_Eigen_DIR_SEARCHED}")
  unset(EIGEN_INCLUDE_DIR CACHE)
endif()
set(_Eigen_DIR_SEARCHED "${Eigen_DIR}" CACHE INTERNAL "Eigen_DIR used by the last Eigen search")

if(Eigen_DIR)
  find_path(EIGEN_INCLUDE_DIR NAMES signature_of_eigen3_matrix_library
            PATHS ${Eigen_DIR} PATH_SUFFIXES eigen3 include/eigen3 NO_DEFAULT_PATH)
else()
  find_path(EIGEN_INCLUDE_DIR NAMES signature_of_eigen3_matrix_library
            PATH_SUFFIXES eigen3 include/eigen3)
endif()

if(EIGEN_INCLUDE_DIR AND EXISTS "${EIGEN_INCLUDE_DIR}/Eigen/Version")
  file(STRINGS "${EIGEN_INCLUDE_DIR}/Eigen/Version" _eigen_version_lines
       REGEX "#define EIGEN_(MAJOR|MINOR|PATCH)_VERSION")
  string(REGEX REPLACE ".*EIGEN_MAJOR_VERSION[ \t]+([0-9]+).*" "\\1" _eigen_major "${_eigen_version_lines}")
  string(REGEX REPLACE ".*EIGEN_MINOR_VERSION[ \t]+([0-9]+).*" "\\1" _eigen_minor "${_eigen_version_lines}")
  string(REGEX REPLACE ".*EIGEN_PATCH_VERSION[ \t]+([0-9]+).*" "\\1" _eigen_patch "${_eigen_version_lines}")
  set(Eigen_VERSION "${_eigen_major}.${_eigen_minor}.${_eigen_patch}")
endif()

# Eigen_VERSION must be a REQUIRED_VARS entry, not only VERSION_VAR: some
# find_package_handle_standard_args implementations report "found" and skip
# the version check when the version variable is empty, which would let a
# pre-5 Eigen (no Eigen/Version, so Eigen_VERSION unset) pass silently.
find_package_handle_standard_args(Eigen
  REQUIRED_VARS EIGEN_INCLUDE_DIR Eigen_VERSION
  VERSION_VAR   Eigen_VERSION)

# A rejected directory must not stick in the cache, or the next configure
# never searches again even after Eigen_DIR is corrected.
if(NOT Eigen_FOUND)
  set(Eigen_REJECTED_DIR "${EIGEN_INCLUDE_DIR}")
  unset(EIGEN_INCLUDE_DIR CACHE)
endif()
mark_as_advanced(EIGEN_INCLUDE_DIR)
