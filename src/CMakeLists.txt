 ### CMakeLists.txt ---
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - RICAM-Linz.
######################################################################
## 
### Commentary: 
## 
######################################################################
## 
### Change log:
## 
######################################################################


SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR})

FOREACH(subdir ${SUBDIRS})
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/CMakeLists.txt")
  ADD_SUBDIRECTORY(${subdir})
ELSE()
  message("Ignoring directory ${subdir}")
endif()

ENDFOREACH()

if(GISMO_BUILD_LIB)

FOREACH(subdir ${SUBDIRS})

install(DIRECTORY   "${CMAKE_CURRENT_SOURCE_DIR}/${subdir}"
        DESTINATION include/${PROJECT_NAME}/
        FILES_MATCHING PATTERN "*.h"
        )

ENDFOREACH()

endif(GISMO_BUILD_LIB)
