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

import os, sys
gismo_path=os.path.join(os.path.dirname(__file__), "../build/lib")
print("G+Smo path:",gismo_path,"(change if needed).")
sys.path.append(gismo_path)

import pygismo as gs
import numpy as np

c = np.array([[0.,1.],[1.,1.]])
print("Coefficient array:",c)

b = gs.nurbs.gsBSpline(0.0,1.0,0,1,c,0,False)
print("B-Spline:", b)

print("Degree:", b.degree(0))

print("Samples:", b.sample(10))

print("Coefficients:\n", b.coefs())

u = np.array([0.5])
print(f"Evaluation of the Bspline on {u[0]}:\n", b.eval(u))


val = np.empty(2)
b.eval_into(u,val)
print(f"Evaluation of the Bspline on {u[0]} with void function:\n", val)

upts = np.linspace(0,1,5)
'''
eval and eval_into take in input a matrix u of size d x N, where each column of u represents one evaluation point
print(f"Evaluation of the Bspline on {u[x] for x in range(len(upts))}:\n", b.eval(upts))

vals = np.empty([2,5])
b.eval_into(upts,vals)
print(f"Evaluation of the Bspline on {u[x] for x in range(len(upts))} with void function:\n", vals)
'''

b.insertKnot(u, 2)
print("Number of coefficients after knot-insertion:\n", b.numCoefs())

b.degreeElevate(1,0)
print("Augmented Degree: ", b.degree(0))


print(f"First derivatives of the Bspline on {u[0]}:\n", b.deriv(u))
derival = np.empty(2)
b.deriv_into(u, derival)
print(f"First derivatives of the Bspline on {u[0]} with void function:\n", derival)


print(f"Second derivatives of the Bspline on {u[0]}:\n", b.deriv2(u))
deriv2val = np.empty(2)
deriv2val = b.deriv2(u)
print(f"Second derivatives of the Bspline on {u[0]} with void function:\n", deriv2val)

#References to matrices can get invalidated
print("--- Updating a referenced array works")
c=b.coefs()
c[0,0]= 2.0
print("Matrix (g+smo):\n", b.coefs())
print("Matrix (ndarray):\n", c)
print("shape:", c.shape)
print("type:", c.dtype)
b.insertKnot(0.5,1) # knot, multiplicity
print("--- Resizing a referenced matrix invalidates the python ndarray")
print("Matrix (g+smo):\n", b.coefs())
print("Matrix (ndarray):\n", c)
print("shape:", c.shape)
print("type:", c.dtype)
