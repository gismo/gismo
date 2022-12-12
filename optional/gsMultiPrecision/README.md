# Multi-precision floating-point extension

G+Smo extension for the [GMP](https://gmplib.org) and [MPFR](https://www.mpfr.org) libraries.

|CMake flags|```-DGISMO_WITH_GMP=ON``` (default ```OFF```) <br> or <br> ```-DGISMO_WITH_MPFR=ON``` (default ```OFF```)|
|--:|---|
|License|[MPL 2.0](https://www.mozilla.org/en-US/MPL/2.0/)|
|OS support|Linux, macOS, Windows (only MinGW)|
|Status|completed|
|Developer|Angelos Mantzaflaris, Matthias Moller|
|Maintainer|M.Moller@tudelft.nl|
|Last checked|10-12-2020|

***

The multi-precision floating-point extension is a wrapper around the
GNU Multiple Precision Arithmetic Library ([GMP](https://gmplib.org))
and the GNU MPFR Library ([MPFR](https://www.mpfr.org)).

From the GMP library, the exact multiprecision rational number type
`mpq_class` can be used as G+Smo's primary data type for coefficients
by setting `-DGISMO_COEFF_TYPE=mpq_class`.

From the MPFR library, the multi-precision real number type
`mpfr::mpreal` can be used as G+Smo's primary data type for
coefficients by setting `-DGISMO_COEFF_TYPE=mpfr::mpreal`.

If `-DGISMO_COEFF_TYPE=mpq_class` or `-DGISMO_COEFF_TYPE=mpfr::mpreal`
is set then `-DGISMO_WITH_GMP=ON` and/or `-DGISMO_WITH_MPFR=ON` are
set automatically.

Note that `-DGISMO_COEFF_TYPE=mpq_class` is incompatible with
`-DGISMO_INDEX_TYPE=int64_t` and `-DGISMO_INDEX_TYPE=long long`.
