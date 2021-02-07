# CoDiPack extension

G+Smo extension for the [CoDiPack - Code Differentiation Package](https://www.scicomp.uni-kl.de/software/codi/).

|CMake flags|```-DGISMO_WITH_CODIPACK=ON``` (default ```OFF```)|
|--:|---|
|Required additional CMake flags|```-DCMAKE_CXX_STANDARD=11``` (or better)<br>```-DGISMO_BUILD_LIB=OFF```|
|License|[MPL 2.0](https://www.mozilla.org/en-US/MPL/2.0/)|
|OS support|Linux, Windows, macOS|
|Status|completed|
|Developer|Matthias Möller|
|Maintainer|M.Moller@tudelft.nl|
|Last checked|10-12-2020|

***
__Table of content__
1. [Introduction](#introduction)
2. [Usage example](#usage_example)
***

__Introdution__

The CoDiPack extension builds on the open-source
[CoDiPack](https://github.com/SciCompKL/CoDiPack) (Code
Differentiation Package) library developed by the [Scientific
Computing Group](http://www.scicomp.uni-kl.de/) at the TU
Kaiserslautern. CoDiPack is a tool for gradient evaluation in computer
programs. It supports the features:

-  Forward mode of Algorithmic Differentiation (AD)
-  Reverse mode of Algorithmic Differentiation (AD)
-  Different tape implementations
-  An AdjointMPI interface
-  External functions
-  Higher order derivatives

G+Smo's CoDiPack extension provides the scalar types
```codi::RealForward```, ```codi::RealReverse```, and
```codi::RealReverseUnchecked``` and several others. The
```codi::RealForward``` type implements the forward mode of AD and the
```codi::RealReverse``` type implements the reverse mode of AD. The
third type, ```codi::RealReverseUnchecked``` is also an implementation
of the reverse mode of AD but it should only be used by experienced
users since it disabled internal checks. For each type there is also a
type with single precession, i.e. ```codi::RealForwardFloat```,
```codi::RealReverseFloat```, and
```codi::RealReverseUncheckedFloat```. The above scalar types can be
used as drop-in replacements of ```real_t``` to enable algorithmic
differentiation (AD) in forward and reverse mode, respectively.

This extension requires C++11 (```-DCMAKE_CXX_STANDARD=11```) or
better enabled and compilation in header-only mode
(```-DGISMO_BUILD_LIB=OFF```).

__Usage example__

The file ```codipack_example.cpp``` illustrates the basic usage of the CoDiPack extension.

1.  Configuration and compilation
    ```bash
    mkdir build
    cd build
    cmake .. -DCMAKE_CXX_STANDARD=11 -DGISMO_BUILD_LIB=OFF -DGISMO_WITH_CODIPACK=ON
    make codipack_example -j4
    ```
2.  Execution
    ```bash
    ./bin/codipack_example
    ```
