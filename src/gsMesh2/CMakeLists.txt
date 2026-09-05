######################################################################
## CMakeLists.txt --- gsMesh2
## This file is part of the G+Smo library.
##
## surface meshes and subdivision
######################################################################

gismo_add_module(gsMesh2
  DEPENDS gsCore gsIO gsUtils gsNurbs gsContainers
  SUBDIRS gsSubdivisionSchemes)
