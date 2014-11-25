 ### CMakeLists.txt ---
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - RICAM-Linz.
######################################################################

## #################################################################
## Collect Headers
## #################################################################


# Collect header file names
aux_header_directory(${CMAKE_CURRENT_SOURCE_DIR} FILES)


set(ALL_HEADERS ${ALL_HEADERS}
    ${FILES}
   CACHE INTERNAL "${PROJECT_NAME} header files"
)
