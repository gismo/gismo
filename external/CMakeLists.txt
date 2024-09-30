######################################################################
### CMakeLists.txt --- external
## This file is part of the G+Smo library.
## 
## Author: Angelos Mantzaflaris 
######################################################################

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

#include_directories(${CMAKE_CURRENT_SOURCE_DIR})
