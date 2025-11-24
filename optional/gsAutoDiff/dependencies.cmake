# dependencies.cmake for gsAutoDiff module
# This file is automatically included by the main CMakeLists.txt
# when gsAutoDiff is in GISMO_OPTIONAL

# Find MUMPS
find_package(autodiff REQUIRED)

if(autodiff_FOUND)
  message(STATUS "AutoDiff found: ${autodiff_INCLUDE_DIR}")
  message(STATUS "AutoDiff libraries: ${autodiff_LIBRARIES}")

  # Add AutoDiff include directories to global includes
  set(GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${autodiff_INCLUDE_DIR}
    CACHE INTERNAL "${PROJECT_NAME} include directories")

  # Add AutoDiff libraries to global linker
  set(gismo_LINKER ${gismo_LINKER} ${autodiff_LIBRARIES}
    CACHE INTERNAL "${PROJECT_NAME} extra linker objects")

else()
  message(WARNING "MUMPS not found. gsAutoDiff module will be disabled.")
endif()