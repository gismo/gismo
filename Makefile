######################################################################
## Makefile
## This file is part of the G+Smo library.
##
## Uses CMake to configure and build G+Smo in the subfolder ./build
##
## Author: Angelos Mantzaflaris
######################################################################

SHELL := /bin/bash
RM    := rm -rf
MKDIR := mkdir -p

ifndef VERBOSE
MAKEFLAGS += --no-print-directory
endif

all: ./build/Makefile
	@ $(MAKE) -C build  && echo "Build completed in folder ./build"

./build/Makefile:
	@ $(MKDIR) build
	@ (cd build >/dev/null 2>&1 && cmake ..)

cleancmake:
	@- $(RM) CMakeCache.txt
	@- $(RM) CMakeFiles
	@- $(RM) cmake_install.cmake

distclean:
	@- (cd build >/dev/null 2>&1 && cmake .. >/dev/null 2>&1)
	@- $(MAKE) --silent -C build clean || true
	@- $(RM) ./build/Makefile
	@- $(RM) ./build/CMake*
	@- $(RM) ./build/cmake.*
	@- $(RM) ./build/*.cmake

ifeq ($(findstring distclean,$(MAKECMDGOALS)),)
ifeq ($(findstring cleancmake,$(MAKECMDGOALS)),)

    $(MAKECMDGOALS): ./build/Makefile
	@ $(MAKE) -C build $(MAKECMDGOALS) && echo "Build completed in folder ./build"

endif
endif
