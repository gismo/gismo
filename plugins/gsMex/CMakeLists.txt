
CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )

PROJECT( gsMex )

##############################################################################
# FIXME: THESE SHOULD BE SET AUTOMATICALLY
##############################################################################
# Specify MEX file generator/program, extension of MEX files, and
# installation directory for MATLAB files.
SET( MEX_GENERATOR /opt/Matlab/bin/mex)
SET( MEX_EXTENSION .mexglx )
SET( MAT_INSTALL_DIR $ENV{HOME}/Matlab )

# Ensure existence of MEX generator.
FIND_PROGRAM( MEX_EXECUTABLE ${MEX_GENERATOR} )
IF( NOT MEX_EXECUTABLE )
  MESSAGE( SEND_ERROR "Failed to find MEX generator: ${MEX_GENERATOR}." )
ENDIF( NOT MEX_EXECUTABLE )

# Find the G+Smo library.
find_package(gismo REQUIRED)

##############################################################################
# FIXME: we should actually include all Gismo include directories
# instead of only the first (element 0 of GISMO_INCLUDE_DIR)
##############################################################################
# Set each directory to be included
list( GET GISMO_INCLUDE_DIRS 0 GISMO_INCLUDE_DIR )
set( GISMO_INCLUDE_EXTERNAL ${gismo_SOURCE_DIR}/external)
set( GISMO_INCLUDE_BINARY ${gismo_BINARY_DIR})
set( GISMO_LIBRARY_DIR L${gismo_BINARY_DIR}/lib )

# Fetch all .cpp files to be mex'ed, and loop over these to set up
# mex-targets.
FILE( GLOB CXX_SRCS src/*.cpp )
SET(MEX_SRCS)
SET(MEX_OUTS)
SET(MEX_TGTS)
FOREACH(f ${CXX_SRCS})
  GET_FILENAME_COMPONENT(fnam ${f} NAME_WE)
  SET(fout ${CMAKE_CURRENT_BINARY_DIR}/${fnam})
  SET(fmex ${CMAKE_CURRENT_BINARY_DIR}/${fnam}${MEX_EXTENSION})
  ADD_CUSTOM_COMMAND(
    OUTPUT  "${fout}"
    COMMAND ${MEX_EXECUTABLE} 
      -cxx
      -I${GISMO_INCLUDE_DIR}
      -I${GISMO_INCLUDE_EXTERNAL}
      -I${GISMO_INCLUDE_BINARY}
      -I${gsMex_SOURCE_DIR}/include
      -L${CMAKE_BINARY_DIR}/lib -lgismo
      -output ${fout}
      ${f}
    DEPENDS ${f} gismo
    COMMENT "Generating ${fmex}"
  )
  LIST(APPEND MEX_SRCS ${fout})
  LIST(APPEND MEX_OUTS ${fmex})
  LIST(APPEND MEX_TGTS ${fnam})
  ADD_CUSTOM_TARGET(${fnam} DEPENDS ${fout})
ENDFOREACH()

ADD_CUSTOM_TARGET(gsMex ALL DEPENDS ${MEX_SRCS})

# Set up installation process.
IF( NOT MAT_INSTALL_DIR )
  MESSAGE( FATAL_ERROR "MATLAB install directory NOT found." )
ENDIF( NOT MAT_INSTALL_DIR )
INSTALL ( FILES ${MEX_OUTS} # MEX files
          DESTINATION "${MAT_INSTALL_DIR}" )
INSTALL ( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/m/ # M files
          DESTINATION "${MAT_INSTALL_DIR}" 
          FILES_MATCHING PATTERN "*.m"
          PATTERN ".svn" EXCLUDE )
