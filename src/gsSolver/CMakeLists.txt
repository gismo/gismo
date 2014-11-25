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

set(gismo_HEADERS ${gismo_HEADERS}
    ${FILES}
   CACHE INTERNAL "${PROJECT_NAME} header files"
)

# Collect template header file names
aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR} HPPFILES)

set(gismo_TMPL_HEADERS ${gismo_TMPL_HEADERS}
    ${HPPFILES}
   CACHE INTERNAL "${PROJECT_NAME} template header files"
)

# Collect source file names
if ( GISMO_BUILD_SHARED_LIB)
   aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} CPPFILES)
   set(gismo_SOURCES ${gismo_SOURCES}
       ${CPPFILES}
       CACHE INTERNAL "${PROJECT_NAME} source files"
)
endif()

## #################################################################
## Generate instance files
## #################################################################

if(GISMO_BUILD_SHARED_LIB)
aux_generate_instance( ${CMAKE_CURRENT_SOURCE_DIR})
endif(GISMO_BUILD_SHARED_LIB)
