#!/usr/bin/python

""""
    @file BSpline curve example

    @brief Play with a B-spline curve in Python

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore
"""

import sys
sys.path.append('/user/amantzaf/home/root/pybind11/build/lib')
 
import pygismo as gs
import numpy as np

c = np.array([[0.,1.],[1.,1.]])
print("Coefficient array:",c)

b = gs.nurbs.gsBSpline(0.0,1.0,0,1,c,0,False)
print("B-Spline:", b)

print("Degree:", b.degree(0))

print("Samples:", b.sample(10))

print("Coefficients:", b.coefs())
