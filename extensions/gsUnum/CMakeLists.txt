######################################################################
### CMakeLists.txt --- gsUnum
## This file is part of the G+Smo library.
## 
## Author: Matthias Moller
######################################################################

## Unum extension
project(gsUnumExtension)

# Apply same configuration as G+Smo
include(gsConfig)

# Fetch Unum sources at configure time
include(gsFetch)
gismo_fetch_directory(Unum
    GIT_REPOSITORY https://github.com/stillwater-sc/universal.git
    DESTINATION    external
)

set(UNUM_INCLUDE_DIR ${gismo_SOURCE_DIR}/external/Unum/posit CACHE INTERNAL "")
include_directories(${UNUM_INCLUDE_DIR})

install(DIRECTORY ${UNUM_INCLUDE_DIR}
        DESTINATION include/posit/ OPTIONAL)
