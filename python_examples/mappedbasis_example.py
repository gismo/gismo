#!/usr/bin/python

""""
    @file mappedbasis_example.py

    @brief Construction and evaluation of a mapped basis

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
"""

import os, sys
gismo_path=os.path.join(os.path.dirname(__file__), "../build/lib")
print("G+Smo path:",gismo_path,"(change if needed).")
sys.path.append(gismo_path)

import pygismo as gs
import numpy as np
import scipy.sparse as sparse

## See gismo/filedata/surfaces/simple.xml for the geometry
c1 = np.array([0.,0.,0.,0.5,1.,1.,1.])
c2 = np.array([0.,0.,0.,0.5,1.,1.,1.])
ku1 = gs.nurbs.gsKnotVector(c1,3)
ku2 = gs.nurbs.gsKnotVector(c2,3)

tbasis = gs.nurbs.gsTensorBSplineBasis2(ku1,ku2)

coefs1 = np.array([
                    [0     ,0    ,0   ],
                    [0.25  ,0    ,0   ],
                    [0.75  ,0    ,0   ],
                    [1     ,0    ,0   ],
                    [0     ,0.25 ,0   ],
                    [0.25  ,0.25 ,0   ],
                    [0.75  ,0.25 ,0   ],
                    [1     ,0.25 ,0   ],
                    [0     ,0.50 ,0   ],
                    [0.25  ,0.50 ,0   ],
                    [0.75  ,0.50 ,0   ],
                    [1     ,0.50 ,0   ],
                    [0     ,0.75 ,0   ],
                    [0.25  ,0.75 ,0   ],
                    [0.75  ,0.75 ,0   ],
                    [1     ,0.75 ,0   ],
                    [0     ,1    ,0   ],
                    [0.25  ,1    ,0   ],
                    [0.75  ,1    ,0   ],
                    [1     ,1    ,0   ],
                        ])

coefs2 = np.array([
                    [1     ,0    ,0   ],
                    [1.25  ,0    ,0   ],
                    [1.75  ,0    ,0   ],
                    [2     ,0    ,0   ],
                    [1     ,0.25 ,0   ],
                    [1.25  ,0.25 ,0   ],
                    [1.75  ,0.25 ,0   ],
                    [2     ,0.25 ,0   ],
                    [1     ,0.50 ,0   ],
                    [1.25  ,0.50 ,0   ],
                    [1.75  ,0.50 ,0   ],
                    [2     ,0.50 ,0   ],
                    [1     ,0.75 ,0   ],
                    [1.25  ,0.75 ,0   ],
                    [1.75  ,0.75 ,0   ],
                    [2     ,0.75 ,0   ],
                    [1     ,1    ,0   ],
                    [1.25  ,1    ,0   ],
                    [1.75  ,1    ,0   ],
                    [2     ,1    ,0   ],
                        ])

tspline1 = gs.nurbs.gsTensorBSpline2(tbasis,coefs1)
tspline2 = gs.nurbs.gsTensorBSpline2(tbasis,coefs2)

mp = gs.core.gsMultiPatch()
mp.addPatch(tspline1)
mp.addPatch(tspline2)

row = np.array([0,   1,   4,   5,   8,   9,   12,  13,  2,   3,   16,  3,   16,  17,  6,   7,   20,  7,   20,  21,  10,  11,  24,  11,  24,  25,  14,  15,  28,  15,  28,  29,  18,  19,  22,  23,  26,  27,  30,  31])
col = np.array([0,   1,   2,   3,   4,   5,   6,   7,   8,   8,   8,   9,   9,   9,   10,  10,  10,  11,  11,  11,  12,  12,  12,  13,  13,  13,  14,  14,  14,  15,  15,  15,  16,  17,  18,  19,  20,  21,  22,  23])
data= np.array([1,   1,   1,   1,   1,   1,   1,   1,   1,   0.5, 0.5, 0.5, 0.5, 1,   1,   0.5, 0.5, 0.5, 0.5, 1,   1,   0.5, 0.5, 0.5, 0.5, 1,   1,   0.5, 0.5, 0.5, 0.5, 1,   1,   1,   1,   1,   1,   1,   1,   1])

mb = gs.core.gsMultiBasis(mp)

mat = sparse.csr_matrix((data, (row, col)), shape=(32, 24))

mbasis = gs.msplines.gsMappedBasis2(mb,mat)