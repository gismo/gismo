## gsVersion.cmake
## 
## Author: Angelos Mantzaflaris 
######################################################################


##########################################
# G+Smo version configuration
#
# Versioning scheme:
#  
#  M.m.b (major.minor.build)
#
#  0.8 will be used for the first alpha release
#      Subsequent updates will be denoted 0.8.1, 0.8.2, ..
#
#  0.9 will be used for the first beta release
#      Subsequent updates will be denoted 0.9.1, 0.9.2, ..
#
#  1.0 will be used for the first milestone (stable) release
##########################################

set(${PROJECT_NAME}_MAJOR_VERSION "0")
set(${PROJECT_NAME}_MINOR_VERSION "8")
set(${PROJECT_NAME}_BUILD_VERSION "0")
set(${PROJECT_NAME}_VERSION_TYPE  "alpha")

# Constuct version info
if( "${PROJECT_NAME}_MINOR_BUILD" STREQUAL "0")
    set(${PROJECT_NAME}_VERSION
        "${${PROJECT_NAME}_MAJOR_VERSION}.${${PROJECT_NAME}_MINOR_VERSION}${${PROJECT_NAME}_VERSION_TYPE}"
        CACHE INTERNAL "${PROJECT_NAME} version number")
else()
    set(${PROJECT_NAME}_VERSION
        "${${PROJECT_NAME}_MAJOR_VERSION}.${${PROJECT_NAME}_MINOR_VERSION}.${${PROJECT_NAME}BUILD_VERSION}${${PROJECT_NAME}_VERSION_TYPE}"
        CACHE INTERNAL "${PROJECT_NAME} version number")
endif()


##########################################
# Try to detect SVN revision
##########################################

IF (EXISTS ${PROJECT_SOURCE_DIR}/.svn) 
    #Determine revision level 
    FIND_PACKAGE(Subversion 1.6) 
    IF(Subversion_FOUND)
       Subversion_WC_INFO(${PROJECT_SOURCE_DIR} GISMO) 
    ELSE(Subversion_FOUND) 
      set(GISMO_WC_REVISION "??")
    ENDIF(Subversion_FOUND)
ELSE (EXISTS ${PROJECT_SOURCE_DIR}/.svn)
    IF (EXISTS ${PROJECT_SOURCE_DIR}/../.svn) 
       #Determine revision level 
          FIND_PACKAGE(Subversion 1.6) 
          IF(Subversion_FOUND)
            Subversion_WC_INFO(${PROJECT_SOURCE_DIR}/../ GISMO) 
          ELSE(Subversion_FOUND) 
            set(GISMO_WC_REVISION "??")
          ENDIF(Subversion_FOUND)
ELSE (EXISTS ${PROJECT_SOURCE_DIR}/../.svn)
     set(GISMO_WC_REVISION "??")
ENDIF(EXISTS ${PROJECT_SOURCE_DIR}/../.svn) 
ENDIF(EXISTS ${PROJECT_SOURCE_DIR}/.svn) 
