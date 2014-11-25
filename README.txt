======================================================================
=====                       Gismo's README                       =====
======================================================================

Gismo is a pure template C++ library for Isogeometric analysis.
Please refer to

./doc/manual/Gismo_manual.pdf 

for detailed documentation. A doxygen documentation is also available.

Dependancies:

- STL library (if you have gcc then you have it already)


Optional dependancies:

- Qt   (if you want to build the Qt application)

- Axel  (if you want to build the Axel plugin)



Licence.

Authors:
Angelos Mantzaflaris
Gabor Kiss


======================================================================
=====                     Compilation options                    ===== 
======================================================================

GISMO_BUILD_SHARED_LIBS ON
GISMO_BUILD_EXAMPLES OFF
GISMO_BUILD_TESTS ON
GISMO_BUILD_QT_APP ON
GISMO_BUILD_DOX OFF
GISMO_BUILD_PDF OFF


======================================================================
=====             Short description of the source tree           ===== 
======================================================================

The sources are found at:

./Gismo/src/

Some tests of the available objects are found in

./tests/

and some more advanced examples of usage in

./examples/

Upon compilation, an instance of the library is compiled, using double
arithmetic. The instantization is done in the folder

./instance/




======================================================================
