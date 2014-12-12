## #################################################################
## Installation
## #################################################################

message ("  CMAKE_INSTALL_PREFIX    ${CMAKE_INSTALL_PREFIX}")

#set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY true)

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR     lib     CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR     bin     CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")

# Set CMake installation directory
if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR ${INSTALL_LIB_DIR}/cmake)
else()
  #set(DEF_INSTALL_CMAKE_DIR ${INSTALL_LIB_DIR}/cmake/gismo)
   set(DEF_INSTALL_CMAKE_DIR ${INSTALL_LIB_DIR})
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH
  "Installation directory for CMake files")

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

# Add all targets to the build-tree export set
if(GISMO_BUILD_SHARED_LIB)
export(TARGETS gismo
  FILE "${PROJECT_BINARY_DIR}/gismoTargets.cmake" APPEND)
endif()

#if(GISMO_WITH_ONURBS)
#  export(TARGETS gsOpennurbsExtension
#    FILE "${PROJECT_BINARY_DIR}/gismoTargets.cmake" APPEND)
#endif(GISMO_WITH_ONURBS)

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE gismo)

# Get relative paths
#file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"
#   "${INSTALL_INCLUDE_DIR}/${PROJECT_NAME}")
#file(RELATIVE_PATH REL_LIB_DIR "${INSTALL_CMAKE_DIR}"
#   "${INSTALL_LIB_DIR}")


# Create the gismoConfig.cmake and gismoConfigVersion.cmake files

# ... for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/src" 
                      "${PROJECT_SOURCE_DIR}/external"
                      "${PROJECT_BINARY_DIR}"
                      "${PROJECT_SOURCE_DIR}/extensions" )
set(CONF_LIB_DIRS     "${CMAKE_BINARY_DIR}/lib")
set(CONF_USE_FILE     "${PROJECT_SOURCE_DIR}/cmake/gismoUse.cmake")
configure_file(${PROJECT_SOURCE_DIR}/cmake/gismoConfig.cmake.in
              "${CMAKE_BINARY_DIR}/gismoConfig.cmake" @ONLY)

# ... for the install tree
set(CONF_INCLUDE_DIRS "${INSTALL_INCLUDE_DIR}/${PROJECT_NAME}")
set(CONF_LIB_DIRS     "${INSTALL_LIB_DIR}")
set(CONF_USE_FILE     "${INSTALL_CMAKE_DIR}/gismoUse.cmake")
configure_file(${PROJECT_SOURCE_DIR}/cmake/gismoConfig.cmake.in
              "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/gismoConfig.cmake" @ONLY)

# ... for both
configure_file(${PROJECT_SOURCE_DIR}/cmake/gismoConfigVersion.cmake.in
  "${CMAKE_BINARY_DIR}/gismoConfigVersion.cmake" @ONLY)

if(GISMO_BUILD_SHARED_LIB)

set_target_properties(gismo PROPERTIES
  PUBLIC_HEADER "${PROJECT_SOURCE_DIR}/Gismo/gismo.h"
  #;${CMAKE_CURRENT_BINARY_DIR}/config.h"
  )

# For gsExport.h
install(FILES ${PROJECT_BINARY_DIR}/gsCore/gsExport.h
        DESTINATION include/${PROJECT_NAME}/gsCore/)

# For gsLinearAlgebra.h
install(DIRECTORY ${PROJECT_SOURCE_DIR}/external/Eigen
        DESTINATION include/${PROJECT_NAME}/ 
        PATTERN "*.txt" EXCLUDE
        PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ
        )

# For gsCmdLine.h
install(DIRECTORY ${PROJECT_SOURCE_DIR}/external/tclap
        DESTINATION include/${PROJECT_NAME}/ 
        FILES_MATCHING
        PATTERN "*.h"
        PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)

# For gsTemplate.h
install(FILES ${PROJECT_SOURCE_DIR}/external/eiquadprog.hpp
        DESTINATION include/${PROJECT_NAME}/ )

# For gsXmlUtils.h
install(FILES ${PROJECT_SOURCE_DIR}/external/rapidxml/rapidxml.hpp
              ${PROJECT_SOURCE_DIR}/external/rapidxml/rapidxml_print.hpp	
        DESTINATION include/${PROJECT_NAME}/rapidxml/)


# For pure install
#install(DIRECTORY ${PROJECT_SOURCE_DIR}/external/rapidxml
#        DESTINATION include/${PROJECT_NAME}/
#        FILES_MATCHING
#        PATTERN "*.hpp"
#        PATTERN ".svn" EXCLUDE
#        PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)

# For gsConfig.h
install(FILES ${PROJECT_BINARY_DIR}/gsCore/gsConfig.h
        DESTINATION include/${PROJECT_NAME}/gsCore/)

# Install cmake files
install(FILES
  "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/gismoConfig.cmake"
  "${CMAKE_BINARY_DIR}/gismoConfigVersion.cmake"
  "${PROJECT_SOURCE_DIR}/cmake/gismoUse.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)
 
# Install the export set for use with the install-tree
#install(EXPORT gismoTargets DESTINATION
#  "${INSTALL_CMAKE_DIR}" COMPONENT dev)

else(GISMO_BUILD_SHARED_LIB)
   message ("Configure with -DGISMO_BUILD_SHARED_LIB=ON to compile the library")
endif(GISMO_BUILD_SHARED_LIB)
