""""
    @file Pythonizations for cppyy bindings

    @brief Improves the usability of the Python bindings

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Scholz
"""

# Pythonizations for the gismo namespace
def pythonize_gismo(klass, name):
    # Enable conversion to numpy
    if name.startswith("gsMatrix<double"):
        def tonumpy(self):
            import numpy as np
            rows = self.rows()
            cols = self.cols()
            numpyarray = np.frombuffer(self.data(), dtype=np.double, count=rows*cols).reshape(rows, cols)
            return numpyarray
        #try:
            setattr(klass, "tonumpy", tonumpy)
        #except KeyError:
        #    pass
