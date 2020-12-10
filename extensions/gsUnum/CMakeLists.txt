######################################################################
### CMakeLists.txt --- gsUnum
## This file is part of the G+Smo library.
## 
## Author: Matthias Moller
######################################################################

# Universal number arithmetic (Unum) extension
project(gsUnumExtension)

if(CMAKE_CXX_STANDARD EQUAL "98" OR CMAKE_CXX_STANDARD EQUAL "11" )
  message(FATAL_ERROR "CMAKE_CXX_STANDARD must be 14 or better for Unum extensions to compile.")
endif()

# Apply same configuration as G+Smo
include(gsConfig)

##
## Universal number arithmetic library
##

# Download Unum sources at configure time
include(gsFetch)
gismo_fetch_directory(unum
    GIT_REPOSITORY https://github.com/stillwater-sc/universal.git
    DESTINATION    external
)

# Set Unum include directory
set(UNUM_INCLUDE_DIR ${gismo_SOURCE_DIR}/external/unum/include CACHE INTERNAL "")

# Add Unum include directory to G+Smo standard include directories
set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${UNUM_INCLUDE_DIR}
  CACHE INTERNAL "gismo include directories")

# Install Unum header files
install(DIRECTORY ${UNUM_INCLUDE_DIR}
        DESTINATION include/gismo/unum/ OPTIONAL)
