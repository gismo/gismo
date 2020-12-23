############################################################################
#
#  Program:         SuperLU
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   October 2, 1995
#
#  Modified:        February 4, 1997  Version 1.0
#                   November 15, 1997 Version 1.1
#                   September 1, 1999 Version 2.0
#
############################################################################
#
#  The machine (platform) identifier to append to the library names
#
PLAT = _@CMAKE_SYSTEM_NAME@

#
#  The name of the libraries to be created/linked to
#
SuperLUroot  = @CMAKE_CURRENT_BINARY_DIR@/superlu
SUPERLULIB   = $(SuperLUroot)/lib/libsuperlu_4.3.a

# Use BLAS library found by CMake
BLASDEF      = -DUSE_VENDOR_BLAS
BLASLIB      = @BLAS_LIBRARIES@ -pthreads

TMGLIB       = libtmglib.a
LIBS         = $(SUPERLULIB) $(BLASLIB)

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         = @CMAKE_AR@
ARCHFLAGS    = cr
RANLIB       = @CMAKE_RANLIB@

CC           = @CMAKE_C_COMPILER@
CFLAGS       = @CMAKE_C_FLAGS_RELEASE@
NOOPTS       = 
FORTRAN      = @CMAKE_Fortran_COMPILER@
FFLAGS       = @CMAKE_Fortran_FLAGS_RELEASE@
LOADER       = $(CC)
LOADOPTS     = @CMAKE_EXE_LINKER_FLAGS_RELEASE@

#
#  C preprocessor defs for compilation for the Fortran interface
#  (-DNoChange, -DAdd_, -DAdd__, or -DUpCase)
#
CDEFS        = -DAdd_
#
# The directory in which Matlab is installed
#
MATLAB       = /usr/sww/matlab

