######################################################################
## gsEigen.cmake
## This file is part of the G+Smo library.
##
## Locates Eigen (>= 5.0, header-only). In order of preference:
##
##  1. Eigen_DIR      a plain Eigen source tree (contains Eigen/Core),
##                    e.g. an unpacked release; used as is
##  2. an installed Eigen3 package, find_package(Eigen3 5.0 CONFIG)
##  3. the release tarball GISMO_EIGEN_VERSION fetched into external/Eigen
##                    (GISMO_FETCH_EIGEN, default ON; the directory is
##                    ignored by git and reused by later configure runs)
##
## Sets GISMO_EIGEN_INCLUDE_DIR (the directory holding Eigen/ and
## unsupported/) and GISMO_EIGEN_SOURCE (Eigen_DIR | installed | fetched).
## gsLibrary.cmake and gsInstall.cmake use these to export the include
## directory or the Eigen3 package dependency to consumers.
######################################################################

set(Eigen_DIR "" CACHE PATH "Path to an Eigen source tree (>= 5.0, containing Eigen/Core). Empty: use an installed Eigen3 package or fetch the release")
option(GISMO_FETCH_EIGEN "Fetch the Eigen release when no installed Eigen3 package is found" ON)
set(GISMO_EIGEN_VERSION "5.0.0" CACHE STRING "Eigen release to fetch when no installed Eigen3 package is found")
set(GISMO_EIGEN_URL "https://gitlab.com/libeigen/eigen/-/archive/${GISMO_EIGEN_VERSION}/eigen-${GISMO_EIGEN_VERSION}.tar.gz"
  CACHE STRING "Download location of the Eigen release")
set(GISMO_EIGEN_SHA256 "315c881e19e17542a7d428c5aa37d113c89b9500d350c433797b730cd449c056"
  CACHE STRING "SHA256 of the Eigen release tarball (empty: not checked)")

if(Eigen_DIR)
  if(NOT EXISTS "${Eigen_DIR}/Eigen/Core")
    message(FATAL_ERROR "Eigen_DIR='${Eigen_DIR}' does not point to an Eigen source tree (expected ${Eigen_DIR}/Eigen/Core).")
  endif()
  set(GISMO_EIGEN_INCLUDE_DIR "${Eigen_DIR}")
  set(GISMO_EIGEN_SOURCE "Eigen_DIR")
else()
  find_package(Eigen3 5.0 CONFIG QUIET)
  if(Eigen3_FOUND)
    get_target_property(GISMO_EIGEN_INCLUDE_DIR Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
    set(GISMO_EIGEN_SOURCE "installed")
  elseif(GISMO_FETCH_EIGEN)
    if(GISMO_EIGEN_SHA256)
      set(_eigen_hash URL_HASH SHA256=${GISMO_EIGEN_SHA256})
    else()
      set(_eigen_hash)
    endif()
    if(NOT COMMAND gismo_fetch_directory)
      include(gsFetch)
    endif()
    gismo_fetch_directory(Eigen
      URL ${GISMO_EIGEN_URL} ${_eigen_hash}
      DESTINATION external)
    if(NOT EXISTS "${gismo_SOURCE_DIR}/external/Eigen/Eigen/Core")
      message(FATAL_ERROR "Fetching Eigen ${GISMO_EIGEN_VERSION} failed. Install Eigen (>= 5.0), or point Eigen_DIR to an Eigen source tree.")
    endif()
    set(GISMO_EIGEN_INCLUDE_DIR "${gismo_SOURCE_DIR}/external/Eigen")
    set(GISMO_EIGEN_SOURCE "fetched")
  else()
    message(FATAL_ERROR "No Eigen (>= 5.0) found. Install Eigen, point Eigen_DIR to an Eigen source tree, or set GISMO_FETCH_EIGEN=ON.")
  endif()
endif()
set(GISMO_EIGEN_INCLUDE_DIR "${GISMO_EIGEN_INCLUDE_DIR}" CACHE INTERNAL "Directory holding Eigen/ and unsupported/")
set(GISMO_EIGEN_SOURCE "${GISMO_EIGEN_SOURCE}" CACHE INTERNAL "Where Eigen comes from: Eigen_DIR, installed or fetched")
message(STATUS "Eigen: ${GISMO_EIGEN_SOURCE} (${GISMO_EIGEN_INCLUDE_DIR})")
