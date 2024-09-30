######################################################################
### CMakeLists.txt --- gsUniversal
## This file is part of the G+Smo library.
## 
## Author: Matthias Moller
######################################################################

# Universal number arithmetic (Unum) extension
project(gsUniversal)
add_custom_target(${PROJECT_NAME})

if(CMAKE_CXX_STANDARD EQUAL "98" OR CMAKE_CXX_STANDARD EQUAL "11" )
  message(FATAL_ERROR "CMAKE_CXX_STANDARD must be 14 or better for Unum extensions to compile.")
endif()

# Apply same configuration as G+Smo
include(gsConfig)

##
## Universal number arithmetic library
##

# Download Unum sources at configure time
#include(gsFetch)
gismo_fetch_directory(universal
    GIT_REPOSITORY https://github.com/stillwater-sc/universal.git
    GIT_TAG main
    DESTINATION    external
)

# Set Unum include directory
set(UNIVERSAL_INCLUDE_DIR ${gismo_SOURCE_DIR}/external/universal/include CACHE INTERNAL "")

# Add Unum include directory to G+Smo standard include directories
set (GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} ${UNIVERSAL_INCLUDE_DIR}
  CACHE INTERNAL "gismo include directories")

# Install Unum header files
install(DIRECTORY ${UNIVERSAL_INCLUDE_DIR}/universal
        DESTINATION include/gismo/)

# Install gsUniversal extension header files
install(DIRECTORY ${PROJECT_SOURCE_DIR}
        DESTINATION include/gismo/
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
