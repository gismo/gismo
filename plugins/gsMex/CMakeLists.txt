######################################################################
## CMakeLists.txt --- gsMex
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
######################################################################

project( gsMex )

#For El Capitan
#see https://gitlab.kitware.com/cmake/community/wikis/doc/cmake/RPATH-handling
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Note: set Matlab_ROOT_DIR for specific installation
# set MATLAB_FIND_DEBUG to debug the find procedure
find_package(Matlab REQUIRED COMPONENTS MAIN_PROGRAM MEX_COMPILER) #MX_LIBRARY
message("Matlab_ROOT_DIR ${Matlab_ROOT_DIR}")
message("Matlab_MEX_COMPILER ${Matlab_MEX_COMPILER}")
message("Matlab_MEX_EXTENSION ${Matlab_MEX_EXTENSION}")

#octave
#set( Matlab_MEX_COMPILER "mkoctfile --mex" CACHE STRING "path to Mex generator")

#############
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
  set(fmex ${CMAKE_CURRENT_BINARY_DIR}/${fnam}.${Matlab_MEX_EXTENSION})

  #matlab_add_mex(-- todo
  add_custom_command(OUTPUT "${fmex}"
      COMMAND ${Matlab_MEX_COMPILER} #-v 
      -cxx
      -I${GISMO_INCLUDE_DIR}
      -I${GISMO_INCLUDE_EXTERNAL}
      -I${GISMO_INCLUDE_BINARY}
      -I${gsMex_SOURCE_DIR}/include
      -L${CMAKE_BINARY_DIR}/lib -lgismo
      -output ${fout} ${f}
      #&& install_name_tool -change @rpath/libgismo.0.dylib ${CMAKE_BINARY_DIR}/lib/libgismo.0.dylib ${fmex}
      DEPENDS ${f} src/mex_common.h
      COMMENT "Building MEX object ${fmex}" )

  list(APPEND MEX_OUTS ${fmex})

endforeach()

add_custom_target(${PROJECT_NAME} ALL DEPENDS gismo ${MEX_OUTS})

set_target_properties(${PROJECT_NAME} PROPERTIES LABELS "${PROJECT_NAME}" FOLDER "${PROJECT_NAME}")

add_custom_command(TARGET ${PROJECT_NAME}
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E echo 'To get started in Matlab execute set_path.m from  ${CMAKE_CURRENT_BINARY_DIR}'
  COMMAND ${CMAKE_COMMAND} -E echo 'Tests are located at ${CMAKE_CURRENT_SOURCE_DIR}/tests'
  )
      
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/set_path.m.in" "${CMAKE_CURRENT_BINARY_DIR}/set_path.m")

INSTALL ( FILES ${MEX_OUTS} # MEX files
          DESTINATION "${Matlab_ROOT_DIR}/bin" )

INSTALL ( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/m/ # M files
          DESTINATION "${Matlab_ROOT_DIR}/bin" 
          FILES_MATCHING PATTERN "*.m")
