
# Gismo version configuration

IF (EXISTS ${PROJECT_SOURCE_DIR}/.svn) 
    #Determine revision level 
    FIND_PACKAGE(Subversion 1.6) 
    IF(Subversion_FOUND)
       Subversion_WC_INFO(${PROJECT_SOURCE_DIR} GISMO) 
    ELSE(Subversion_FOUND) 
      set(GISMO_WC_REVISION svn)
    ENDIF(Subversion_FOUND)
ELSE (EXISTS ${PROJECT_SOURCE_DIR}/.svn)
    IF (EXISTS ${PROJECT_SOURCE_DIR}/../.svn) 
       #Determine revision level 
          FIND_PACKAGE(Subversion 1.6) 
          IF(Subversion_FOUND)
            Subversion_WC_INFO(${PROJECT_SOURCE_DIR}/../ GISMO) 
          ELSE(Subversion_FOUND) 
            set(GISMO_WC_REVISION svn)
          ENDIF(Subversion_FOUND)
ELSE (EXISTS ${PROJECT_SOURCE_DIR}/../.svn)
    set(GISMO_WC_REVISION 0)
ENDIF(EXISTS ${PROJECT_SOURCE_DIR}/../.svn) 
ENDIF(EXISTS ${PROJECT_SOURCE_DIR}/.svn) 

# Set version data
set(${PROJECT_NAME}_MAJOR_VERSION "0" CACHE STRING "gismo major version number.")
set(${PROJECT_NAME}_MINOR_VERSION "0" CACHE STRING "gismo minor version number.")
set(${PROJECT_NAME}_VERSION
  "${${PROJECT_NAME}_MAJOR_VERSION}.${${PROJECT_NAME}_MINOR_VERSION}.${GISMO_WC_REVISION}")
mark_as_advanced(${PROJECT_NAME}_MAJOR_VERSION)
mark_as_advanced(${PROJECT_NAME}_MINOR_VERSION)

add_definitions(-DGISMO_VERSION="${${PROJECT_NAME}_VERSION}")
