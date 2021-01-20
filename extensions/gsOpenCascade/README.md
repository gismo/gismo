# Open CASCADE Technology extension

G+Smo extension for the [Open CASCADE Technology](https://dev.opencascade.org) (OCCT) software development kit.

|||
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
packages to be installed on your system. These are essentially OpenGL,
TCL/TK and X11. Below we give instructions for different operating
systems that are known to work.

__Linux__

_CentOS/Red Hat_

1.  Installation of the general development tools
    ```bash
    sudo yum group install "Development Tools"
    ```
2.  Installation of the additional libraries and header files
    ```bash
    sudo yum install freetype-devel libXi-devel libXmu-devel mesa-libGL-devel tk-devel
    ```

_Debian/Ubuntu_

1.  Installation of the general development tools
    ```bash
    sudo apt-get install build-essential
    ```
2.  Installation of the additional libraries and header files
    ```bash
    sudo apt-get install libgl-dev libxi-dev libxmu-dev mesa-common-dev tk-dev
    ```

__macOS__

_Windows_
