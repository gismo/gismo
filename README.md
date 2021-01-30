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
=====                     version 20.12 Alpha                    =====
=====                   https://github.com/gismo                 =====
======================================================================
```

# Continuous Integration status
| **System** | **Status** | **More information** |
|------------|------------|----------------------|
| [CDash](https://cdash-ci.inria.fr/index.php?project=Gismo) | [![cdash](https://img.shields.io/website?down_color=lightgrey&down_message=offline&label=CDash&up_color=green&up_message=up&url=https%3A%2F%2Fcdash-ci.inria.fr%2Findex.php%3Fproject%3DGismo)](https://cdash-ci.inria.fr/index.php?project=Gismo) | Report results from all builds |
| [Appveyor](https://ci.appveyor.com/project/gismo/gismo)  | [![Appveyor status](https://ci.appveyor.com/api/projects/status/abps59xbt1gjwci1/branch/stable?svg=true)](https://cdash-ci.inria.fr/index.php?project=Gismo&filtercount=1&field1=site&compare1=63&value1=[appVeyor]) | Windows MSVC 14.0 |
| [Circle CI](https://circleci.com/gh/gismo/gismo) | [![Circle CI](https://circleci.com/gh/gismo/gismo.svg?style=svg)](https://cdash-ci.inria.fr/index.php?project=Gismo&filtercount=1&field1=site&compare1=63&value1=[cci]) | macOS XCode9, C++98 <br> macOS XCode10, C++11 <br> macOS XCode11, C++14 <br> macOS XCode12, C++17 |
| [Codeship](https://app.codeship.com/projects/123289)  | [![Codeship Status](https://app.codeship.com/projects/2aa19360-8998-0133-39fd-66416d65b267/status?branch=stable)](https://cdash-ci.inria.fr/index.php?project=Gismo&filtercount=1&field1=site&compare1=63&value1=[codeship]) | |
| [GitLab](https://gitlab.com/gismo-ci/gismo/-/pipelines)    | [![pipeline status](https://gitlab.com/gismo-ci/gismo/badges/gitlab_ci/pipeline.svg)](https://cdash-ci.inria.fr/index.php?project=Gismo&filtercount=1&field1=site&compare1=63&value1=[gitlab-ci]) | Linux GCC6, C++98 <br> Linux GCC7, C++11 <br> Linux GCC8, C++14 <br> Linux GCC9, C++17 <br> Linux GCC10, C++20 <br> Linux Clang7, C++98 <br> Linux Clang8, C++11 <br> Linux Clang9, C++14 <br> Linux Clang10, C++17 <br> Linux Clang11, C++20] |
| [Travis](https://travis-ci.org/gismo/gismo/branches) | [![Travis Status](https://travis-ci.org/gismo/gismo.svg?branch=stable)](https://cdash-ci.inria.fr/index.php?project=Gismo&filtercount=1&field1=site&compare1=63&value1=[travis]) | macOS XCode9 C++98 <br> macOS XCode10, C++11 <br> macOS XCode11, C++14 <br> Linux GCC, C++98 |
| [GitHub Actions](https://github.com/gismo/gismo/actions) | [![Build Status](https://github.com/gismo/gismo/workflows/gismo/badge.svg?branch=stable)](https://cdash-ci.inria.fr/index.php?project=Gismo&filtercount=1&field1=site&compare1=63&value1=[actions]) | Latest Linux/MacOS/Windows |
| [Jenkins](https://ci.inria.fr/gismo/job/gismo/job/gismo/job/stable) | [![Build Status](https://ci.inria.fr/gismo/buildStatus/icon?job=gismo%2Fgismo%2Fstable)](https://cdash-ci.inria.fr/index.php?project=Gismo&filtercount=1&field1=site&compare1=63&value1=[jenkins]) |VMs for Linux/MacOS/Windows |
| GCC Farm | [Status](https://cdash-ci.inria.fr/index.php?project=Gismo&filtercount=1&field1=site&compare1=63&value1=[gccfarm]) | Builders from the GCC Farm   |
| [OBS](https://build.opensuse.org/package/show/home:filiatra/gismo) | [binaries](https://software.opensuse.org/download/package?project=home:filiatra&package=gismo)  | Upstream package builds for many Linux distributions |
| [Launchpad](https://code.launchpad.net/~g+smo/+recipe/g+smo-daily) |[binaries](https://launchpad.net/~g+smo/+archive/ubuntu/upstream/+packages)  | Upstream package builds for Ubuntu distributions |


This README file contains brief information. More details are found in
the [Wiki pages](https://github.com/gismo/gismo/wiki).

The latest revision of the code can be obtained using git (via https):

```git clone https://github.com/gismo/gismo.git```

or using subversion:

```svn co https://github.com/gismo/gismo/trunk gismo```

or as a zip file:

https://github.com/gismo/gismo/archive/stable.zip

# Prerequisites

* Operating systems:
  - MS Windows
  - Linux
  - macOS

* Configuration: [CMake 2.8.12](https://cmake.org) or newer.

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

More information at https://github.com/gismo/gismo/wiki

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

  https://github.com/gismo/gismo/wiki

* Bug reports:

  https://github.com/gismo/gismo/issues


# People

Coordinator and maintainer: Angelos Mantzaflaris

See full list in [our wiki pages](https://github.com/gismo/gismo/wiki/About--G-Smo)

# OS-license

The G+Smo library is distributed under the Mozilla Public License v2.0.  (see [LICENSE.txt](https://github.com/gismo/gismo/blob/stable/LICENSE.txt)).

