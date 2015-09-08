######################################################################
## CMakeLists.txt ---
## This file is part of the G+Smo library. 
##
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - 2015 RICAM-Linz.
######################################################################

#include (GenerateExportHeader)

#include(cotire)
#cotire(${${PROJECT_NAME}_MODULES})

## #################################################################
## Add library targets
## #################################################################

#message("Using ${${PROJECT_NAME}_EXTENSIONS}")
#message("Using ${${PROJECT_NAME}_MODULES}")
#message("Using ${${PROJECT_NAME}_SOURCES}")

if(GISMO_BUILD_LIB)

#if ("x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xGNU")
#  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-implicit-templates")
#endif()

if("x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xMSVC")
 set(${PROJECT_NAME}_SOURCES ${${PROJECT_NAME}_SOURCES} 
     "${gismo_SOURCE_DIR}/src/misc/gsDllMain.cpp")
endif()

  add_library(${PROJECT_NAME} SHARED
    ${${PROJECT_NAME}_MODULES}
    ${${PROJECT_NAME}_SOURCES}
    ${${PROJECT_NAME}_EXTENSIONS}
    )

  set_target_properties(${PROJECT_NAME} PROPERTIES 
  PUBLIC_HEADER "${PROJECT_SOURCE_DIR}/src/${PROJECT_NAME}.h" 
  POSITION_INDEPENDENT_CODE ON
  #COMPILE_DEFINITIONS ${PROJECT_NAME}_EXPORTS # Used for DLL exporting (defined by default by CMake)
  FOLDER "G+Smo libraries"
  )

  if (GISMO_WITH_ONURBS AND "x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xMSVC" )
    target_link_libraries(${PROJECT_NAME} Rpcrt4)
  endif()

  if (GISMO_WITH_PSOLID)
    target_link_libraries(${PROJECT_NAME} ${PARASOLID_LIBRARY})
  endif()

  if (GISMO_WITH_SUPERLU)
    target_link_libraries(${PROJECT_NAME} ${SUPERLU_LIBRARIES})
  endif()

if (GISMO_WITH_PARDISO)
   if (PARDISO_USE_MKL)
     find_package(MKL REQUIRED)
     target_link_libraries(${PROJECT_NAME} ${MKL_LIBRARIES})
   else()
     find_package(Pardiso REQUIRED)
     target_link_libraries(${PROJECT_NAME} Pardiso)
   endif()
endif()

#  if (GISMO_WITH_OPENMP)
#    find_package(OpenMP REQUIRED)
#  endif()

  if (GISMO_WITH_IPOPT)
     #might be empty before download
     target_link_libraries(${PROJECT_NAME} ${IPOPT_LIBRARY})
  endif()

  if (GISMO_EXTRA_DEBUG AND DBGHELP_FOUND) 
    target_link_libraries(${PROJECT_NAME} ${DBGHELP_LIBRARY}) 	
  endif() 

 set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)

 if( WIN32 ) # Copy the dll to the bin folder to allow executables to find it
    if(CMAKE_CONFIGURATION_TYPES)
      add_custom_command(
      TARGET ${PROJECT_NAME}
      POST_BUILD
      #COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/bin
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/bin/$<CONFIGURATION>
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:${PROJECT_NAME}> ${CMAKE_BINARY_DIR}/bin/$<CONFIGURATION>
      COMMAND ${CMAKE_COMMAND} -E echo 'The file $<TARGET_FILE:${PROJECT_NAME}> is copied to the bin folder for convenience.' )
    else()
      add_custom_command(
      TARGET ${PROJECT_NAME}
      POST_BUILD
      #COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/bin
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/bin
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:${PROJECT_NAME}> ${CMAKE_BINARY_DIR}/bin
      COMMAND ${CMAKE_COMMAND} -E echo 'The file $<TARGET_FILE:${PROJECT_NAME}> is copied to the bin folder for convenience.' )
    endif()
  endif( WIN32 )

endif(GISMO_BUILD_LIB)

  add_library(${PROJECT_NAME}_static STATIC
  ${${PROJECT_NAME}_MODULES}
  ${${PROJECT_NAME}_SOURCES}
  ${${PROJECT_NAME}_EXTENSIONS}
  )

  set_target_properties(${PROJECT_NAME}_static PROPERTIES 
  POSITION_INDEPENDENT_CODE ON
  EXCLUDE_FROM_ALL 1
  EXCLUDE_FROM_DEFAULT_BUILD 1
  COMPILE_DEFINITIONS ${PROJECT_NAME}_STATIC
  LABELS ${PROJECT_NAME}
  )
 
  if (MKL_FOUND)
    # See http://eigen.tuxfamily.org/dox/TopicUsingIntelMKL.html
    target_link_libraries(${PROJECT_NAME} ${MKL_LIBRARIES})
  endif(MKL_FOUND)

  IF (GISMO_EXTRA_DEBUG AND DBGHELP_FOUND) 
     TARGET_LINK_LIBRARIES(${PROJECT_NAME}_static ${DBGHELP_LIBRARY}) 	
  ENDIF() 

  set_target_properties(${PROJECT_NAME}_static 
  PROPERTIES LINKER_LANGUAGE CXX
  OUTPUT_NAME ${PROJECT_NAME}_static )


set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib/)


## #################################################################
## Installation
## #################################################################

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR     lib     CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR     bin     CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")


if(GISMO_BUILD_LIB)

  install(TARGETS ${PROJECT_NAME}
  # IMPORTANT: Add the ${PROJECT_NAME} library to the "export-set"
  EXPORT gismoTargets
  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/${PROJECT_NAME}"
  #ARCHIVE DESTINATION lib
  )

endif(GISMO_BUILD_LIB)
