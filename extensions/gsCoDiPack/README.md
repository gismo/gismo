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

The CoDiPack extension provides two new scalar types
`codi::RealForward` and `codi::RealReverse` that can be used in place
of `real_t` to enable algorithmic differentiation (AD) in forward and
reverse mode, respectively. This extension requires C++11 or better
and compilation in header-only mode.
