######################################################################
## CMakeLists.txt --- gsOptimizer
## This file is part of the G+Smo library.
##
## optimization problems and solvers
######################################################################

if(GISMO_WITH_IPOPT)
  include_directories(${IPOPT_INCLUDE_DIR})
endif()

gismo_add_module(gsOptimizer
  DEPENDS gsCore gsIO gsUtils)
