######################################################################
## dependencies.cmake --- gsSpectra
## This file is part of the G+Smo library.
##
## Included by the top-level CMakeLists.txt before add_subdirectory(src):
## core includes <gsSpectra/gsSpectra.h> under gsSpectra_ENABLED, so the
## Spectra include directory must be in GISMO_INCLUDE_DIRS before any core
## target is created.
######################################################################

gismo_add_dependency(Spectra
  FIND_PACKAGE   spectra
  GIT_REPOSITORY https://github.com/yixuan/spectra.git
  GIT_TAG        v1.2.0
  MODE           HEADER_ONLY
  INCLUDE_SUBDIR include
  TARGET         Spectra::Spectra)
