
project( gsMex )

#find_package(Matlab)

# Specify MEX file generator/program, extension of MEX files, and
# installation directory for MATLAB files.
SET( MEX_GENERATOR mex CACHE STRING "path to Mex generator")

if ( CMAKE_SYSTEM_NAME MATCHES "Linux" )
  if (CMAKE_SIZEOF_VOID_P EQUAL 8)
    SET( MEX_EXTENSION mexa64 CACHE STRING "Mex extension")
  else()
    SET( MEX_EXTENSION mexa64 CACHE STRING "Mex extension")
  endif()
elseif ( CMAKE_SYSTEM_NAME MATCHES "Windows" )
  SET( MEX_EXTENSION mex264 CACHE STRING "Mex extension")
elseif ( CMAKE_SYSTEM_NAME MATCHES "Darwin" )
  SET( MEX_EXTENSION mexmaci64 CACHE STRING "Mex extension")
else ()
  SET( MEX_EXTENSION mex-UNKNOWN CACHE STRING "Mex extension")
endif()

# Ensure existence of MEX generator.
FIND_PROGRAM( MEX_EXECUTABLE ${MEX_GENERATOR} )
IF( NOT MEX_EXECUTABLE )
  MESSAGE( SEND_ERROR "Failed to find MEX generator: ${MEX_GENERATOR}." )
ENDIF( NOT MEX_EXECUTABLE )

##############################################################################
# FIXME: we should actually include all Gismo include directories
# instead of only the first (element 0 of GISMO_INCLUDE_DIR)
# foreach(...)
##############################################################################
list( GET GISMO_INCLUDE_DIRS 0 GISMO_INCLUDE_DIR )
set( GISMO_INCLUDE_EXTERNAL ${gismo_SOURCE_DIR}/external)
set( GISMO_INCLUDE_BINARY ${gismo_BINARY_DIR})

# Fetch all .cpp files to be mex'ed, and loop over these to set up
# mex-targets.
file( GLOB CXX_SRCS src/*.cpp)
set(MEX_OUTS)
foreach(f ${CXX_SRCS})
  GET_FILENAME_COMPONENT(fnam ${f} NAME_WE)
  set(fout ${CMAKE_CURRENT_BINARY_DIR}/${fnam})
  set(fmex ${CMAKE_CURRENT_BINARY_DIR}/${fnam}.${MEX_EXTENSION})

  add_custom_command(OUTPUT "${fmex}"
      COMMAND ${MEX_EXECUTABLE} 
      -cxx
      -I${GISMO_INCLUDE_DIR}
      -I${GISMO_INCLUDE_EXTERNAL}
      -I${GISMO_INCLUDE_BINARY}
      -I${gsMex_SOURCE_DIR}/include
      -L${CMAKE_BINARY_DIR}/lib -lgismo
      -output ${fout} ${f}
      #&& install_name_tool -change @rpath/libgismo.0.dylib ${CMAKE_BINARY_DIR}/lib/libgismo.0.dylib ${fmex}
      DEPENDS ${f} include/mex_common.h
      COMMENT "Building MEX object ${fmex}" )

  list(APPEND MEX_OUTS ${fmex})

endforeach()

add_custom_target(${PROJECT_NAME} ALL DEPENDS gismo ${MEX_OUTS})
set_target_properties(${PROJECT_NAME} PROPERTIES LABELS "${PROJECT_NAME}" FOLDER "${PROJECT_NAME}")

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/set_path.m.in" "${CMAKE_CURRENT_BINARY_DIR}/set_path.m")

# Set up installation process.
get_filename_component(MAT_INSTALL_DIR ${MEX_GENERATOR} DIRECTORY)
IF( NOT MAT_INSTALL_DIR )
  MESSAGE( FATAL_ERROR "MATLAB install directory NOT found." )
ENDIF( NOT MAT_INSTALL_DIR )

INSTALL ( FILES ${MEX_OUTS} # MEX files
          DESTINATION "${MAT_INSTALL_DIR}" )

INSTALL ( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/m/ # M files
          DESTINATION "${MAT_INSTALL_DIR}" 
          FILES_MATCHING PATTERN "*.m")
