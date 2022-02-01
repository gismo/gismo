#!/usr/bin/python

""""
    @file fileIO_example.py

    @brief File input/output using pygismo

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
"""

import os, sys
gismo_path=os.path.join(os.path.dirname(__file__), "../build/lib")
print("G+Smo path:",gismo_path,"(change if needed).")
sys.path.append(gismo_path)

import pygismo as gs
import numpy as np
import scipy.sparse as sparse


"""" First example: save data for bvp """

# [!Geometry]
c1 = np.array([0.,0.,1.,1.])
c2 = np.array([0.,0.,1.,1.])
ku1 = gs.nurbs.gsKnotVector(c1,1)
ku2 = gs.nurbs.gsKnotVector(c2,1)

coefs = np.array([
                    [0     ,0  ],
                    [1     ,0  ],
                    [0     ,1  ],
                    [1     ,1  ],
                        ])


# Construct basis using knot vectors
tbasis1 = gs.nurbs.gsTensorBSplineBasis2(ku1,ku2)
tspline1 = gs.nurbs.gsTensorBSpline2(tbasis1,coefs)

print("Coefficients:\n", tspline1.coefs())

mp = gs.core.gsMultiPatch()
mp.addPatch(tspline1)
mp.addPatch(tspline1)

""" Other way: read geometry from given file """
#mp = gs.core.gsMultiPatch()
#file = gs.io.gsFileData("planar/square_with_disk.xml")
#file.getAnyFirst(mp) # Assume that there exist only one gsMultiPatch
# [!Geometry]

# [!Right hand side]
f = gs.core.gsFunctionExpr("256*pi*pi*pi*pi*(4*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2)
# [!Right hand side]

# [!Exact solution]
ms = gs.core.gsFunctionExpr("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",2)
# [!Exact solution]

# [!Boundary]
dirichlet = gs.core.gsFunctionExpr("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",2)
neumann = gs.core.gsFunctionExpr(" -4*pi*(cos(4*pi*y) - 1)*sin(4*pi*x)",
                                 " -4*pi*(cos(4*pi*x) - 1)*sin(4*pi*y)",2)

bcs = gs.pde.gsBoundaryConditions()

#               patch_nr, side, boundary condition, function
bcs_info = [[0, gs.core.side.west, gs.pde.bctype.dirichlet, dirichlet],
            [0, gs.core.side.east, gs.pde.bctype.dirichlet, dirichlet],
            [0, gs.core.side.south, gs.pde.bctype.dirichlet, dirichlet],
            [0, gs.core.side.north, gs.pde.bctype.dirichlet, dirichlet],
            [0, gs.core.side.west, gs.pde.bctype.neumann, neumann],
            [0, gs.core.side.east, gs.pde.bctype.neumann, neumann],
            [0, gs.core.side.south, gs.pde.bctype.neumann, neumann],
            [0, gs.core.side.north, gs.pde.bctype.neumann, neumann]]

for bc_pair in bcs_info:
    #               patch_nr, side, boundary condition, function, unknown, parametric, component
    bcs.addCondition(bc_pair[0],bc_pair[1],bc_pair[2],bc_pair[3],0,False,0)

#bcs.setGeoMap(mp)
# [!Boundary]

# [!Option list]
opt = gs.io.gsOptionList()
opt.addSwitch("plot", "Plotting the results.", True)
# [!Option list]

# [!Save the data to the XML-file]
file = gs.io.gsFileData()
file.add(bcs)   # id=0 Boundary
file.add(f)     # id=1 Source function
file.add(opt)   # id=2 Optionlist
file.add(ms)    # id=3 Exact solution
file.add(mp)    # id=4+#patches Geometry (should be last!)
file.save("test_bvp.xml", False)
print("Filedata saved: test_bvp.xml")

## In G+Smo:
    # gsFileData<real_t> fd("test_bvp.xml");
    # gsMultiPatch<> mp_file;
    # gsBoundaryConditions<> bcInfo_file;
    # gsFunctionExpr<> sourceFunction_file;
    # gsOptionList opList_file;
    # gsFunctionExpr<> exactSol_file;
    #
    # fd.getId(0, bcInfo_file);
    # fd.getId(1, sourceFunction_file);
    # fd.getId(2, opList_file);
    # fd.getId(3, exactSol_file);
    # fd.getAnyFirst(mp_file);
    #
    # gsInfo << "mp_file: " << mp_file << "\n";
    # gsInfo << "bcInfo_file: " << bcInfo_file << "\n";
    # gsInfo << "sourceFunction_file: " << sourceFunction_file << "\n";
    # gsInfo << "opList_file: " << opList_file << "\n";
    # gsInfo << "exactSol_file: " << exactSol_file << "\n";
# [!Save the data to the XML-file]

"""" Second example: save multibasis and sparsematrix """
# [!MultiBasis]
mb = gs.core.gsMultiBasis(mp)
# [!MultiBasis]

# [!SparseMatrix]
row = np.array([0,   1,   4,   5,   8,   9,   12,  13,  2,   3,   16,  3,   16,  17,  6,   7,   20,  7,   20,  21,  10,  11,  24,  11,  24,  25,  14,  15,  28,  15,  28,  29,  18,  19,  22,  23,  26,  27,  30,  31])
col = np.array([0,   1,   2,   3,   4,   5,   6,   7,   8,   8,   8,   9,   9,   9,   10,  10,  10,  11,  11,  11,  12,  12,  12,  13,  13,  13,  14,  14,  14,  15,  15,  15,  16,  17,  18,  19,  20,  21,  22,  23])
data= np.array([1,   1,   1,   1,   1,   1,   1,   1,   1,   0.5, 0.5, 0.5, 0.5, 1,   1,   0.5, 0.5, 0.5, 0.5, 1,   1,   0.5, 0.5, 0.5, 0.5, 1,   1,   0.5, 0.5, 0.5, 0.5, 1,   1,   1,   1,   1,   1,   1,   1,   1])

mat = sparse.csr_matrix((data, (row, col)), shape=(32, 24))
# [!SparseMatrix]

# [!Save the data to the XML-file]
file2 = gs.io.gsFileData()
file2.add(mb)
file2.add(mat)
file2.save("test_mspline.xml", False)
print("Filedata saved: test_mspline.xml")
# [!Save the data to the XML-file]



