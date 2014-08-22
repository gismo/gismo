##############################################################
#    CMake Project Wrapper Makefile               
#
# Builds G+Smo in ./build/
#
############################################################## 

SHELL := /bin/bash
RM    := rm -rf
MKDIR := mkdir

all: ./build/Makefile
	@ $(MAKE) -C build  && echo "Build completed in folder ./build"

./build/Makefile:
	@ $(MKDIR) build
	@ (cd build >/dev/null 2>&1 && cmake ..)

distclean:
	@- (cd build >/dev/null 2>&1 && cmake .. >/dev/null 2>&1)
	@- $(MAKE) --silent -C build clean || true
	@- $(RM) ./build/Makefile
	@- $(RM) ./build/CMake*
	@- $(RM) ./build/cmake.*
	@- $(RM) ./build/*.cmake

ifeq ($(findstring distclean,$(MAKECMDGOALS)),)

    $(MAKECMDGOALS): ./build/Makefile
	@ $(MAKE) -C build $(MAKECMDGOALS) && echo "Build completed in folder ./build"

endif
