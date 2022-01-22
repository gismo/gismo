# Unum extension

G+Smo extension for the [Universal Number Arithmetic](https://github.com/stillwater-sc/universal) library.

|CMake flags|```-DGISMO_WITH_UNUM=ON``` (default ```OFF```)|
|--:|---|
|Required additional CMake flags|```-DCMAKE_CXX_STANDARD=14``` (or better)|
|License|[MPL 2.0](https://www.mozilla.org/en-US/MPL/2.0/)|
|OS support|Linux, Windows, macOS|
|Status|completed|
|Developer|Matthias Moller|
|Maintainer|M.Moller@tudelft.nl|
|Last checked|10-12-2020|

***
__Table of content__
1. [Introduction](#introduction)
2. [Posits](#posits)
3. [Usage example](#usage-example)
***

## Introduction

The Universal Number Arithmetic extension builds on the open-source header-only [Universal](https://github.com/stillwater-sc/universal) library developed by Stillwater Supercomputing, Inc. The Universal library provides several drop-in replacements for IEEE floating-point which aim at reproducibility of arithmetic computations especially in parallel execution.

## Posits

G+Smo's Unum extension currently supports the following Posit number systems:

-  `posit_2_0`
-  `posit_3_0`
-  `posit_3_1`
-  `posit_4_0`
-  `posit_8_0`
-  `posit_8_1`
-  `posit_16_1`
-  `posit_32_2`
-  `posit_64_3`
-  `posit_128_4`
-  `posit_256_5`

The above data types can be used as drop-in replacements for `double` and `float` or globally as the G+Smo coefficient type
```bash
cmake .. -DGISMO_COEFF_TYPE=posit_32_2
```

## Usage example
