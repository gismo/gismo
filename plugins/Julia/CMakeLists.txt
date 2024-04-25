######################################################################
## CMakeLists.txt --- Gismo.jl
## This file is part of the G+Smo library.
##
######################################################################

project( Gismo.jl )

#to do
#find_package(Julia)

# Fetch all .cpp files to be mex'ed, and loop over these to set up
# mex-targets.
file( GLOB JL_SRCS src/*.jl)

add_custom_target(${PROJECT_NAME} ALL DEPENDS gismo ${JL_OUTS})
      
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/Gismo_path.jl.in" "${CMAKE_CURRENT_SOURCE_DIR}/src/Gismo_path.jl")

#INSTALL ...

