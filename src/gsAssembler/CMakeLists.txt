######################################################################
## CMakeLists.txt --- gsAssembler
## This file is part of the G+Smo library.
##
## quadrature, expression assembler, visitor assemblers, projections
######################################################################

gismo_add_module(gsAssembler
  DEPENDS gsCore gsMatrix gsUtils gsDomain gsIO gsFunctionExpr gsContainers gsNurbs gsHSplines gsMSplines gsPde gsExpressions gsCInterface)
