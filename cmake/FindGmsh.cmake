# FindGmsh.cmake
# Locates the Gmsh C++ API.
#
# Sets:
#   GMSH_FOUND      - True if Gmsh was found
#   GMSH_INC        - Include directory containing gmsh.h
#   GMSH_LIB        - Path to libgmsh
#   GMSH_INCLUDE_DIRS / GMSH_LIBRARIES - standard aliases
#
# Honors hints (in order):
#   GMSH_INC / GMSH_LIB (cmake vars or env vars)
#   GMSH_DIR (root install prefix) / gmsh_DIR
#   CMAKE_PREFIX_PATH

if(DEFINED ENV{GMSH_INC} AND NOT GMSH_INC)
  set(GMSH_INC $ENV{GMSH_INC})
endif()
if(DEFINED ENV{GMSH_LIB} AND NOT GMSH_LIB)
  set(GMSH_LIB $ENV{GMSH_LIB})
endif()

find_path(GMSH_INC
  NAMES gmsh.h
  HINTS ${GMSH_DIR} ${gmsh_DIR} ENV GMSH_DIR
  PATH_SUFFIXES include
)

find_library(GMSH_LIB
  NAMES gmsh
  HINTS ${GMSH_DIR} ${gmsh_DIR} ENV GMSH_DIR
  PATH_SUFFIXES lib lib64
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Gmsh
  REQUIRED_VARS GMSH_INC GMSH_LIB
)

if(GMSH_FOUND)
  set(GMSH_INCLUDE_DIRS ${GMSH_INC})
  set(GMSH_LIBRARIES ${GMSH_LIB})
endif()

mark_as_advanced(GMSH_INC GMSH_LIB)
