# dependencies.cmake for gsAutoDiff module
#
# WHY THIS FILE EXISTS:
# This file is included by the top-level CMakeLists.txt *before*
# add_subdirectory() is called for any module. That early inclusion is required
# for two reasons:
#
#   1. C++17 standard: autodiff requires C++17. CMAKE_CXX_STANDARD must be
#      raised before gsCore and other targets are created; doing so inside the
#      module's own CMakeLists.txt would be too late.
#
#   2. Global scope: autodiff_FOUND and GISMO_INCLUDE_DIRS are set here so
#      they are visible to all subsequently configured targets, not just to
#      the gsAutoDiff subdirectory.
#
# This file is automatically included by the main CMakeLists.txt
# when gsAutoDiff is in GISMO_OPTIONAL

# The autodiff library (v1.1.2) is C++17-only. G+Smo defaults to C++14
# (see cmake/AddCXXCompileOptions.cmake), which fails to compile the autodiff
# headers. This file is included (from the top-level CMakeLists.txt) after
# gsConfig but before add_subdirectory(src), so raising the standard here
# propagates to gsCore and all subsequently created targets.
if(NOT DEFINED CMAKE_CXX_STANDARD OR CMAKE_CXX_STANDARD LESS 17)
  message(STATUS "gsAutoDiff: autodiff requires C++17; raising CMAKE_CXX_STANDARD to 17")
  set(CMAKE_CXX_STANDARD 17 CACHE INTERNAL "" FORCE)
endif()

# Try system-installed autodiff first
find_package(autodiff QUIET)

if(NOT autodiff_FOUND)
  message(STATUS "AutoDiff not found on system. Fetching v1.1.2 from GitHub...")

  gismo_fetch_directory(autodiff
    GIT_REPOSITORY https://github.com/autodiff/autodiff.git
    GIT_TAG        v1.1.2
    DESTINATION    external
  )

  set(AUTODIFF_FETCH_DIR "${gismo_SOURCE_DIR}/external/autodiff")

  if(EXISTS "${AUTODIFF_FETCH_DIR}/autodiff")
    if(NOT TARGET autodiff::autodiff)
      add_library(autodiff::autodiff INTERFACE IMPORTED GLOBAL)
      target_include_directories(autodiff::autodiff INTERFACE "${AUTODIFF_FETCH_DIR}")
    endif()
    set(autodiff_FOUND TRUE)
    set(autodiff_INCLUDE_DIR "${AUTODIFF_FETCH_DIR}")
    message(STATUS "AutoDiff fetched: ${AUTODIFF_FETCH_DIR}")
  else()
    message(WARNING "Failed to fetch AutoDiff. gsAutoDiff module will be disabled.")
  endif()
endif()

if(autodiff_FOUND)
  message(STATUS "AutoDiff found: ${autodiff_INCLUDE_DIR}")

  set(GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${autodiff_INCLUDE_DIR}
    CACHE INTERNAL "${PROJECT_NAME} include directories")

  set(gismo_LINKER ${gismo_LINKER} ${autodiff_LIBRARIES}
    CACHE INTERNAL "${PROJECT_NAME} extra linker objects")
else()
  message(WARNING "AutoDiff not found. gsAutoDiff module will be disabled.")
endif()
