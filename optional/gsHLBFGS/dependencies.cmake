######################################################################
## dependencies.cmake --- gsHLBFGS
## This file is part of the G+Smo library.
##
## Included by the top-level CMakeLists.txt before add_subdirectory(src):
## core includes <gsHLBFGS/gsHLBFGS.h> (src/gsCore/gsFunction.hpp,
## src/gsModeling/gsBarrierCore.h, gsBarrierPatch.h,
## gsSurfaceReparameterization.h), so the HLBFGS objects must be in
## gismo_EXTENSIONS before the gismo libraries are created.
##
## SOURCE_DIR is the directory CONTAINING HLBFGS/, because gsHLBFGS.h
## includes "HLBFGS/HLBFGS.h". EXPORT_SYMBOLS keeps HLBFGS()/INIT_HLBFGS
## visible in libgismo.so: the header-inline gsHLBFGS<T> calls them from
## consumers of the shared library.
######################################################################

gismo_add_dependency(HLBFGS
  SOURCE_DIR     "${gismo_SOURCE_DIR}/external"
  MODE           SOURCES
  SOURCE_GLOBS   HLBFGS/HLBFGS.cpp HLBFGS/HLBFGS_BLAS.cpp HLBFGS/ICFS.cpp HLBFGS/LineSearch.cpp
  EXPORT_SYMBOLS
  TARGET         HLBFGS::HLBFGS)
