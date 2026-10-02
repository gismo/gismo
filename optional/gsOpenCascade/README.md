# Open CASCADE Technology extension

G+Smo extension for the [Open CASCADE Technology](https://dev.opencascade.org) (OCCT) software development kit.

|GISMO_OPTIONAL| gsOpenCascade|
|--:|---|
|License|[MPL 2.0](https://www.mozilla.org/en-US/MPL/2.0/)|
|OS support|Linux, Windows, macOS|
|Status|completed|
|Developer|Angelos Mantzaflaris|
|Maintainer|angelos.mantzaflaris@inria.fr|
|Last checked|10-12-2020|

***
__Table of content__
1. [Prerequisites](#prerequisites)

***

## Prerequisites

Building the Open CASCADE Technology extension requires additional
packages to be installed on your system. These are essentially Freetype,
OpenGL, TCL/TK and X11 (for visualization support). Below we give
instructions for different operating systems that are known to work.

The CMake option `OCC_ENABLE_X11` (default: `ON`) controls whether the
OpenCascade Visualization module (including `TKV3d` and `TKService`) is
built. When `ON`, the required X11 development headers must be present or
CMake will emit a fatal error with install instructions. Set
`-DOCC_ENABLE_X11=OFF` to build without visualization support.

> **Note for users of the gsGmsh module**: If you also enable `gsGmsh`
> with a system-installed `libgmsh`, that library was compiled against the
> system OpenCASCADE (OCC 7.6.x on Ubuntu 24.04). In that case you **must**
> also install the system OCC development packages so gsOpenCascade uses the
> same ABI rather than downloading an older incompatible version:
> ```bash
> sudo apt-get install libocct-foundation-dev libocct-modeling-algorithms-dev \
>              libocct-modeling-data-dev libocct-data-exchange-dev \
>              libocct-visualization-dev libocct-ocaf-dev
> ```

__Linux__

_CentOS/Red Hat_

1.  Installation of the general development tools
    ```bash
    sudo yum group install "Development Tools"
    ```
2.  Installation of the additional libraries and header files
    ```bash
    sudo yum install freetype-devel libX11-devel libXext-devel libXmu-devel libXi-devel mesa-libGL-devel mesa-libGLU-devel tk-devel
    ```

_Debian/Ubuntu_

1.  Installation of the general development tools
    ```bash
    sudo apt-get install build-essential
    ```
2.  Installation of the additional libraries and header files
    ```bash
    sudo apt-get install libfreetype-dev libx11-dev libxext-dev libxmu-dev libxi-dev libgl-dev mesa-common-dev libglu1-mesa-dev tk-dev
    ```

__macOS__

__Windows__
