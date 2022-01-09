#!/usr/bin/python

""""
    @file biharmonic_example.py

    @brief Compute biharmonic2_example using pygismo

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
"""

import os, sys

gismo_path = os.path.join(os.path.dirname(__file__), "../build/lib")
print("G+Smo path:", gismo_path, "(change if needed).")
sys.path.append(gismo_path)

import pygismo as gs
import numpy as np


# [!Geometry]
mp = gs.core.gsMultiPatch()
file = gs.io.gsFileData("planar/two_squares.xml")
file.getAnyFirst(mp)  # Assume that there exist only one gsMultiPatch

print(mp.nPatches())
boxSide = gs.core.boxSide(gs.core.side.west)
print(boxSide.index()) # get the side index
patchSide = gs.core.patchSide(1,boxSide)
print(patchSide.side().index()) # get the side index
print(patchSide.patch()) # get the patch index
print(mp.boundaries())
for bdy in mp.boundaries():
    print("Patch:", bdy.patch(), "Side:", bdy.side().index())
    print()

# [!Geometry]

# [!Right hand side]
f = gs.core.gsFunctionExpr("256*pi*pi*pi*pi*(4*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))", 2)
# [!Right hand side]

# [!Exact solution]
ms = gs.core.gsFunctionExpr("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)", 2)
# [!Exact solution]

# [!Boundary]
dirichlet = gs.core.gsFunctionExpr("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)", 2)
neumann = gs.core.gsFunctionExpr(" -4*pi*(cos(4*pi*y) - 1)*sin(4*pi*x)",
                                 " -4*pi*(cos(4*pi*x) - 1)*sin(4*pi*y)", 2)

bcs = gs.pde.gsBoundaryConditions()
for bdy in mp.boundaries():
    #               patch_nr, side, boundary condition, function, unknown, parametric, component
    bcs.addCondition(bdy, gs.pde.bctype.dirichlet, dirichlet, 0, False, 0)
    bcs.addCondition(bdy, gs.pde.bctype.neumann, neumann, 0, False, 0)

# bcs.setGeoMap(mp)
# [!Boundary]

# [!Option list]
opt = gs.io.gsOptionList()
opt.addSwitch("plot", "Plotting the results.", True)
# [!Option list]

# [!Save the data to the XML-file]
file = gs.io.gsFileData()
file.add(bcs)  # id=0 Boundary
file.add(f)  # id=1 Source function
file.add(opt)  # id=2 Optionlist
file.add(ms)  # id=3 Exact solution
file.add(mp)  # id=4+#patches Geometry (should be last!)
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
