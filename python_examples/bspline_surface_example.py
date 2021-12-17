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

## See gismo/filedata/surfaces/simple.xml for the geometry
c1 = np.array([0.,0.,0.,0.,1.,1.,1.,1.])
c2 = np.array([0.,0.,0.,0.,1.,1.,1.,1.])
ku1 = gs.nurbs.gsKnotVector(c1,3)
ku2 = gs.nurbs.gsKnotVector(c2,3)

coefs = np.array([
                    [0     ,0    ,0   ],
                    [0.33  ,0    ,0.2 ],
                    [0.66  ,0    ,0.3 ],
                    [1     ,0    ,0.1 ],
                    [0     ,0.33 ,0   ],
                    [0.33  ,0.33 ,0   ],
                    [0.66  ,0.33 ,0.4 ],
                    [1     ,0.33 ,0.2 ],
                    [0     ,0.66 ,0   ],
                    [0.33  ,0.66 ,0.1 ],
                    [0.66  ,0.66 ,0.2 ],
                    [1     ,0.66 ,0.2 ],
                    [0     ,1    ,0.3 ],
                    [0.33  ,1    ,0.2 ],
                    [0.66  ,1    ,0   ],
                    [1     ,1    ,0   ],
                        ])


# Construct basis using knot vectors
tbasis1 = gs.nurbs.gsTensorBSplineBasis2(ku1,ku2)
tspline1 = gs.nurbs.gsTensorBSpline2(tbasis1,coefs)

print("Coefficients:\n", tspline1.coefs())

u = np.array([0.5,0.5])
print(f"Evaluation of the Bspline on {u[0],u[1]}:\n", tspline1.eval(u))


val = np.empty(3)
tspline1.eval_into(u,val)
print(f"Evaluation of the Bspline on {u[0],u[1]} with void function:\n", val)

## NOTE: The lines below give a segmentation fault (probably something wrong with deconstructor?)
# # Construct basis using knot vectors
# bbasis1 = gs.nurbs.gsBSplineBasis(ku1)
# bbasis2 = gs.nurbs.gsBSplineBasis(ku2)
# tbasis2 = gs.nurbs.gsTensorBSplineBasis2(bbasis1,bbasis2)
# tspline2 = gs.nurbs.gsTensorBSpline2(tbasis2,coefs)

# print("Coefficients:\n", tspline2.coefs())

# u = np.array([0.5,0.5])
# print(f"Evaluation of the Bspline on {u[0],u[1]}:\n", tspline2.eval(u))


# val = np.empty(3)
# tspline2.eval_into(u,val)
# print(f"Evaluation of the Bspline on {u[0],u[1]} with void function:\n", val)