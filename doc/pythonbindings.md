Using G+Smo in Python {#python}
======================

G+Smo currently has two separate sets of Python bindings: Bindings based
on [pybind11](https://github.com/pybind/pybind11) and bindings based
on [cppyy](https://cppyy.readthedocs.io/en/latest/). The main difference is 
that the pybind11-bindings are precompiled
while the cppyy-bindings are based on JIT and are generated at runtime.

cppyy bindings
--
Cppyy-bindings are generated at runtime based on just-in-time compilation.
For using G+Smo this is useful, since the bilinear forms are usually assembled
using the expression assembler. Since gismo::gsExprAssembler is based
on templated expressions, its assembly method can only be compiled
once the bilinear form has been defined, which makes precompiled bindings impossible
unless we restrict the available bilinear forms. 
Cppyy's just-in-time compilation solves this problem.

### Configuration
Configure G+Smo with the options `-DGISMO_WITH_CPPYY=On -DCMAKE_CXX_STANDARD=20`.
This creates the Python package `gismo_cppyy` in the subdirectory `cppyy` of your build folder.

Then, build the `gismo` library as usual. 

Run `make wheel` in order to generate the pip wheel in `cppyy/dist`.
It can be installed with `make install-bindings` or directly using `pip`.

The package looks for the G+Smo sources in the source directory and for the shared library in the build directory.
This means that these directories need to stay in place.

### Usage
The bindings can be imported in Python using
```
from gismo_cppyy import gismo
```
All classes and functions from G+Smo can be found in the `gismo` namespace. 

Examples can be found in `python_examples/poisson_example_cppyy.py` and `python_examples/fitting_example_cppyy.py`.

### Numpy compatibility
Data can be exchanged with numpy by converting numpy arrays into gismo objects and vice versa using the
`tonumpy()` method of gismo::gsVector, gismo::gsMatrix and gismo::gsAsMatrix, 
and the respective `fromnumpy()` class methods.

For example, you can run
```
>>> from gismo_cppyy import gismo
>>> import numpy as np
>>> a = np.array([[1,0],[0,1]], dtype=np.double)
>>> A = gismo.gsMatrix["double"].fromnumpy(a)
>>> b = A.tonumpy()
>>> B = gismo.gsAsMatrix["double"].fromnumpy(b)
```

Note that `gismo.gsMatrix["double"].fromnumpy()`
copies the underlying data while `gismo.gsAsMatrix["double"].fromnumpy()` does not.
