######################################################################
## CMakeLists.txt --- gsUtils
## This file is part of the G+Smo library.
##
## utilities (also the legacy gsMesh in gsUtils/gsMesh)
######################################################################

gismo_add_module(gsUtils
  DEPENDS gsParallel
  SUBDIRS gsMesh)
