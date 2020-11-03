```
     GGGGGGGGG      GGGG      GGGGGGGGG  GGGGGG   GGGGGG  GGGGGGGGGG
    GGGG            GGGG     GGGG        GGGGGG  GGGGGG  GGGG   GGGG
   GGGG         GGGGGGGGGGGG GGGGGGGGG   G GGGG  G GGGG GGGG    GGGG
   GGGG GGGGGG GGGGGGGGGGGGG GGGGGGGGGG GG GGGG GG GGGG GGGG   GGGGG
  GGGGG  GGGGG GGGGGGGGGGGG  GGGGGGGGG  GG GGGGGG GGGG  GGGG   GGGG
  GGGG   GGGG      GGGG           GGGG  GG  GGGG  GGGG  GGGG   GGGG
   GGGGGGGGGG      GGGG     GGGGGGGGG  GG   GGG   GGGG  GGGGGGGGGG

======================================================================
=====             Geometry plus Simulation modules               =====
=====                      version 0.8.5 Alpha                   =====
=====                   https://github.com/gismo                 =====
=====          Copyright (C) JKU-RICAM-Linz, 2012 - 2018         =====
======================================================================
```

# Continuous Integration status
| **System** | **Status** | **More information** |
|------------|------------|----------------------|
| CDash      | [![cdash](https://img.shields.io/website?down_color=lightgrey&down_message=offline&label=CDash&up_color=green&up_message=up&url=https%3A%2F%2Fcdash-ci.inria.fr%2Findex.php%3Fproject%3DGismo)](https://cdash-ci.inria.fr/index.php?project=Gismo) | Various builds and configurations on https://cdash-ci.inria.fr/index.php?project=Gismo |
| Linux      | [![Circle CI](https://circleci.com/gh/gismo/gismo.svg?style=svg)](https://circleci.com/gh/gismo/gismo) |
| Linux     | [ ![Codeship Status](https://app.codeship.com/projects/2aa19360-8998-0133-39fd-66416d65b267/status?branch=stable)](https://app.codeship.com/projects/123289)
| macOS      | [![Travis Status](https://travis-ci.org/gismo/gismo.svg?branch=stable)](https://travis-ci.org/gismo/gismo/branches) | 
| Windows    | [![Appveyor status](https://ci.appveyor.com/api/projects/status/abps59xbt1gjwci1/branch/stable?svg=true)](https://ci.appveyor.com/project/gismo/gismo) |

This README file contains brief information. More details are found in
the [Wiki pages](http://gs.jku.at/gismo).

The latest revision of the code can be obtained using subversion:

```svn co https://github.com/gismo/gismo/trunk gismo```

or using git (via https):

```git clone https://github.com/gismo/gismo.git```

or as a zip file:

https://github.com/gismo/gismo/archive/stable.zip

# Prerequisites

* Operating systems:
  - MS Windows
  - Linux
  - macOS

* Configuration: [CMake 2.8.8](https://cmake.org) or newer.

* Compilers tested include recent versions of
  - [AMD Optimizing C/C++ Compiler](https://developer.amd.com/amd-aocc/)
  - [Clang](https://clang.llvm.org) also Apple Clang
  - [GNU GCC](https://gcc.gnu.org)
  - [Intel C++ compiler](https://software.intel.com/content/www/us/en/develop/tools/compilers/c-compilers.html)
  - [Mingw64](http://mingw-w64.org/)
  - [MS Visual Studio C++](https://visualstudio.microsoft.com)
  - [PGI C/C++](https://www.pgroup.com/index.htm) only with `GISMO_WITH_OPENMP=OFF`
  
* Compilers known to not work
  - [Oracle Developer Studio](https://www.oracle.com/application-development/technologies/developerstudio.html) fails to compile Eigen
  - [IBM XLC C/C++](https://www.ibm.com/products/xl-cpp-linux-compiler-power) fails to compile Eigen

* Recommended:
   - [Paraview](https://www.paraview.org) for visualization.

# Compilation

The compilation requires configuration using CMake at a new, empty
folder (in-source builds are disabled).

* On Linux/macOS: A Unix makefile exists in the root source
  folder. Running "make" creates a sub folder named "build" and
  executes CMake and compilation inside that folder. Alternatively,
  choose your own build folder and execute CMake pointing to the
  sources.

* On MS Windows: Run cmake-gui tool (from an environment that is
  configured with your compiler) to generate makefiles (or Visual
  Studio project). Then execute the make tool to launch
  compilation. Alternatively, use the QtCreator GUI and open the
  CMakeLists.txt file on the root folder to create a QtCreator
  project.

After successful compilation a dynamic library is created in ./lib and
executable example programs are output at the ./bin subdirectory of
the build folder.

Additionally, if Doxygen is available on the system one can execute
(eg. on Linux):

```make doc```

to obtain the Doxygen documentation in HTML format. The main doxygen
page is at ./doc/html/index.html.

More information at
http://gs.jku.at/gismo

# Configuration Options

The available options are displayed at CMake configuration.  Short
description and default setting follows:

* CMAKE_BUILD_TYPE        *Release*

  Available values are the standard CMake build configurations: Debug,
Release, RelWithDebInfo, MinSizeRel.

* GISMO_COEFF_TYPE        *double*

  The arithmetic type to be used for all computations. Available options
include double, long double, float.

* GISMO_EXTRA_DEBUG       *OFF*

  If set to ON additional debugging tools are enabled during
compilation. These include checked iterators for GCC and MSVC
compilers and call stack back-trace printout when a runtime exception
occurs.

* GISMO_BUILD_LIB         *ON*

  If enabled a dynamic library is created using GISMO_COEFF_TYPE
arithmetic. A target for a static library named gismo_static is also
created but not compiled by default.

* GISMO_BUILD_EXAMPLES    *ON*

  If enabled the programs in the examples folder are compiled, and
executables are created in build-folder/bin.

* GISMO_BUILD_UNITTESTS   *OFF*

  If enabled the tests in the unittests folder are compiled, and an
executable is created in build-folder/bin.

* GISMO_BUILD_AXL         *OFF*

  If enabled the plugin for Axel modeler is compiled (requires Axel).

* GISMO_WITH_PSOLID       *OFF*

  If enabled the extensions using functionalities of Parasolid geometric
kernel are compiled (requires Parasolid).

* GISMO_WITH_ONURBS       *OFF*

  If enabled the extension for reading and writing of Rhinoceros' 3DM is
compiled.

* CMAKE_INSTALL_PREFIX   (system dependent)

  The location for installation of the library, e.g. /usr/local on some
Linux systems.


# Directory structure


The source tree consists of the following sub-folders:

* **src**

Contains all source files. Code is partitioned into modules. Currently
eleven modules are present as sub-folders:

   - **gsCore**
   - **gsMatrix**
   - **gsNurbs**
   - **gsHSplines**
   - **gsModeling**
   - **gsAssembler**
   - **gsSolver**
   - **gsPde**
   - **gsTensor**
   - **gsIO**
   - **gsUtils**

* **examples**

  Examples of usage, small programs and tutorials.

* **unittests**

  Unittests for some parts of the codebase.

* **filedata**

  Data files in the XML format the G+Smo can read and write.

* **extensions**

  Optional additional features that can be compiled along G+Smo.

* **plugins**

   The plugins for:

   - Axel modeler
   - Rhinoceros' 3DM

* **cmake**

  Cmake configuration files.

* **doc**

  Files related to doxygen documentation.

# Contact and support

* Wiki pages:

  https://gs.jku.at/gismo

* Mailing list: gismo@ricam.oeaw.ac.at

* Subscribe to the mailing list at:

  https://list.ricam.oeaw.ac.at/lists/listinfo/gismo

* Bug reports:

  https://github.com/gismo/gismo/issues


# People

Coordinator and maintainer: Angelos Mantzaflaris

See full list in [our wiki pages](https://github.com/gismo/gismo/wiki/About--G-Smo)

# Acknowledgement and os-license

G+Smo library for isogeometric analysis is supported by the Research
Network ''Geometry + Simulation'' (NFN S117), funded by the Austrian
Science Fund, www.gs.jku.at.  It is distributed under the Mozilla
Public License v2.0.  (see [LICENSE.txt](https://github.com/gismo/gismo/blob/stable/LICENSE.txt)).

