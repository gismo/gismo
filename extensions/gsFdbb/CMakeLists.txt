### CMakeLists.txt ---
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2016 - RICAM-Linz.
######################################################################

## FDBB extension
project(gsFdbbExtension)

# Apply same configuration as G+Smo
include(gsConfig)

gismo_fetch_directory(FDBB
    GIT_REPOSITORY https://gitlab.com/mmoelle1/FDBB.git
    DESTINATION    external
)

set(FDBB_INCLUDE_DIR ${gismo_SOURCE_DIR}/external/FDBB/fdbb CACHE INTERNAL "")
include_directories(${FDBB_INCLUDE_DIR})

install(DIRECTORY ${FDBB_INCLUDE_DIR}
        DESTINATION include/fdbb/ OPTIONAL)
