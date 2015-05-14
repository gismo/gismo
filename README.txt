======================================================================
=====                                                            ===== 
=====                          G+Smo                             =====
=====             Geometry plus Simulation modules               =====
=====                     version 0.8alpha                       =====
=====          Copyright (C) JKU-RICAM-Linz, 2012 - 2015         =====
=====                                                            ===== 
=====                      version 0.8alpha                      ===== 
=====                                                            ===== 
======================================================================

 
G+Smo library for isogeometric analysis is supported by the Research
Network ``Geometry + Simulation'' (NFN S117), funded by the Austrian
Science Fund, www.gs.jku.at.  It is distributed under the Mozilla
Public License v2.0.  (see LICENCE.txt).

This README file contains brief information. More details are found in
the Wiki page: 

https://ricamsvn.ricam.oeaw.ac.at/trac/gismo/wiki

The latest revision of the code can be obtained via SVN (subversion) at:

https://ricamsvn.ricam.oeaw.ac.at/trac/gismo/stable

Username/password: anonsvn

======================================================================
=====                      Prerequisites                         ===== 
======================================================================

* Operating systems: 
  - MS Windows
  - Linux
  - MacOSX

* Configuration: CMake 2.8.8 or newer.

* Compilers tested include recent versions of 
  - GNU GCC
  - MS Visual C++
  - Clang
  - Intel C++ compiler
  - Mingw32 

* Recommended:
   - Paraview for visualization.

======================================================================
=====                        Compilation                         ===== 
======================================================================

The compilation requires configuration using CMake at a new, empty
folder (in-source builds are disabled).

* On Linux/MacOSX: A Unix makefile exists in the root source
  folder. Running "make" creates a sub folder named "build" and
  executes CMake and compilation inside that folder. Alternatively,
  choose your own build folder and execute CMake pointing to the
  sources.

* On MS Windows: Run cmake-gui tool (from an environment that is
  configured with your compiler) o generate makefiles (or Visual
  Studio project). Then execute the make tool to launch
  compilation. Alternatively, use the QtCreator GUI and open the
  CMakeLists.txt file on the root folder to create a QtCreator
  project.

After successful compilation a dynamic library is created in ./lib and
executable example programs are output at the ./bin subdirectory of
the build folder.

Additionally, if Doxygen is available on the system one can execute
(eg. on Linux):

make doc

to obtain the Doxygen documentation in HTML format. The main doxygen
page is at ./doc/html/index.html.

More information at
https://ricamsvn.ricam.oeaw.ac.at/trac/gismo

======================================================================
=====                   Configuration Options                    ===== 
======================================================================

The available options are displayed at CMake configuration.  Short
description and default setting follows:

* CMAKE_BUILD_TYPE        RelWithDebInfo

Available values are the standard CMake build configurations: Debug,
Release, RelWithDebInfo, MinSizeRel.

* GISMO_COEFF_TYPE        double

The arithmetic type to be used for all computations. Available options
include double, long double, float.

* GISMO_EXTRA_DEBUG       OFF

If set to ON additional debugging tools are enabled during
compilation. These include checked iterators for GCC and MSVC
compilers and call stack back-trace printout when a runtime exception
occurs.

* GISMO_BUILD_LIB         ON

If enabled a dynamic library is created using GISMO_COEFF_TYPE
arithmetic. A target for a static library named gismo_static is also
created but not compiled by default.

* GISMO_BUILD_EXAMPLES    ON

If enabled the programs in the examples folder are compiled, and
executables are created in build-folder/bin.

* GISMO_BUILD_AXL         OFF

If enabled the plugin for Axel modeler is compiled (requires Axel).

* GISMO_WITH_PSOLID       OFF

If enabled the extensions using functionalities of Parasolid geometric
kernel are compiled.(requires Parasolid).

* GISMO_WITH_ONURBS       OFF

If enabled the extension for reading and writing of Rhinoceros' 3DM is
compiled.

* CMAKE_INSTALL_PREFIX   (system dependent)

The location for installation of the library, e.g. /usr/local on some
Linux systems.


======================================================================
=====                   Directory structure                      ===== 
======================================================================

The source tree consists of the following sub-folders:

* src

Contains all source files. Code is partitioned into modules. Currently
eleven modules are present as sub-folders:
       
   - gsCore
   - gsMatrix
   - gsNurbs
   - gsHSplines
   - gsModeling
   - gsAssembler
   - gsSolver
   - gsPde
   - gsTensor
   - gsIO
   - gsUtils

* examples

Examples of usage, small programs and tutorials.

* filedata

Data files in the XML format the G+Smo can read and write.

* extensions

Optional additional features that can be compiled along G+Smo.

* plugins

Plugins for exporting functionality of G+Smo to third party software.

* cmake

Cmake configuration files

* doc

Files related to doxygen documentation.

======================================================================
=====                    Contact and support                     ===== 
======================================================================

* Wiki pages:

  https://ricamsvn.ricam.oeaw.ac.at/trac/gismo/wiki

* Mailing list: gismo@ricam.oeaw.ac.at

* Subscribe to the mailing list at: 

  https://list.ricam.oeaw.ac.at/lists/listinfo/gismo

* Bug reports: 

  https://ricamsvn.ricam.oeaw.ac.at/trac/gismo/newticket


======================================================================
=====                             People                         ===== 
======================================================================

* Coordinator and maintainer: Angelos Mantzaflaris

* Developers:

  Andrea Bressan
  Florian Buchegger
  Antonella Falini
  Carlotta Giannelli
  Clemens Hofreither
  Mario Kapl
  Gabor Kiss
  Stefan Kleiss
  Angelos Mantzaflaris
  David Mayer
  Dominik Mokris
  Jarle Sogn
  Jaka Speh
  Stephen Moore
  Manh Dang Nguyen
  Michael Pauley
  Harald Weiner

* Scientific board:

  Bert Juettler
  Ulrich Langer
  Walter Zulehner
