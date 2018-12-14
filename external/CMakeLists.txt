### CMakeLists.txt ---
##
## Author: Angelos Mantzaflaris
## Copyright (C) 2018 - RICAM-Linz.
######################################################################

## #################################################################
## Collect Headers
## #################################################################

set(${PROJECT_NAME}_externals "${CMAKE_CURRENT_SOURCE_DIR}"
    CACHE INTERNAL "Externals directory")

find_package(ZLIB QUIET)
if (NOT ZLIB_FOUND OR GISMO_ZLIB_STATIC)
  set(GISMO_ZLIB_STATIC TRUE CACHE INTERNAL "Use static zlib.")
  add_definitions(-DZ_PREFIX) #use prefixed zlib for gzstream
  add_subdirectory(zlib)
else()
  message(STATUS "Zlib found ${ZLIB_LIBRARIES}")
  set(gismo_LINKER ${gismo_LINKER} ${ZLIB_LIBRARIES}
  CACHE INTERNAL "${PROJECT_NAME} extra linker objects")
  include_directories(SYSTEM ${ZLIB_INCLUDES})
endif()

add_subdirectory(gzstream)

if(GISMO_BUILD_CPPLOT)
#add_subdirectory(boost)
add_subdirectory(cpplot)
endif(GISMO_BUILD_CPPLOT)

#add_subdirectory(Eigen)
#add_subdirectory(rapidxml)

if(GISMO_WITH_MPFR)
  find_package(MPFR QUIET)
  if (NOT MPFR_FOUND  OR GISMO_BUILD_MPFR)
    include(mpfr.cmake)
  else()
    message(STATUS "MPFR found ${MPFR_LIBRARY}")
    set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${MPFR_INCLUDE_DIR}
      CACHE INTERNAL "${PROJECT_NAME} include directories")
    set(gismo_LINKER ${gismo_LINKER} ${MPFR_LIBRARY}
      CACHE INTERNAL "${PROJECT_NAME} extra linker objects")
  endif()
endif()

if(GISMO_WITH_GMP)
  find_package(GMP QUIET)
  if (NOT GMP_FOUND  OR GISMO_BUILD_GMP)
    include(gmp.cmake)
  else()
    message(STATUS "GMP found ${GMP_LIBRARY} ${GMPXX_LIBRARY}")
    set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${GMP_INCLUDE_DIR}
      CACHE INTERNAL "${PROJECT_NAME} include directories")
    set(gismo_LINKER ${gismo_LINKER} ${GMP_LIBRARY} ${GMPXX_LIBRARY}
      CACHE INTERNAL "${PROJECT_NAME} extra linker objects")
  endif()
endif()
