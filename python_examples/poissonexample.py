""""
    @file Poisson equation example

    @brief Solve the Poisson equation using the expression assembler in Python

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Scholz
"""

import math
from gismo import gismo


numElevate = 0
numRefine = 0

fd = gismo.gsFileData["double"]("pde/poisson2d_bvp.xml")
mp = gismo.gsMultiPatch['double']()
#
fd.getId(0, mp)


f = gismo.gsFunctionExpr['double']()
fd.getId(1, f)
print(f"Source function {f}")

bc = gismo.gsBoundaryConditions["double"]()

fd.getId(2, bc)



print("Got bc from file")
bc.setGeoMap(mp)
print(f"Boundary conditions:\n{bc}")

ms = gismo.gsFunctionExpr["double"]()
fd.getId(3, ms)

Aopt = gismo.gsOptionList()

fd.getId(4, Aopt)

dbasis = gismo.gsMultiBasis["double"](mp, True)
dbasis.setDegree(dbasis.maxCwiseDegree() + numElevate)

print(f"Patches: {mp.nPatches()}, degree: {dbasis.minCwiseDegree()}")

A = gismo.gsExprAssembler["double"]()

A.setOptions(Aopt)

print(f"Active Options:\n{A.options()}")

A.setIntegrationElements(dbasis)
ev = gismo.gsExprEvaluator["double"](A)

G = A.getMap(mp)

u = A.getSpace(dbasis)

ff = A.getCoeff(f, G)

u_ex = ev.getVariable(ms, G)


solver = gismo.gsSparseSolver["double"].CGDiagonal()


for r in range(3):
    print(f"ref {r}")
    dbasis.uniformRefine()


    basis = u.source().basis(0)

    u.setup(bc, gismo.dirichlet.l2Projection, 0)

    A.initSystem()


    print(f"Number of DoFs: {A.numDofs()}")

    A.assemble(gismo.expr.igrad(u, G) * gismo.expr.igrad(u,G).tr() * gismo.expr.meas(G), u * ff * gismo.expr.meas(G))
    

    g_N = A.getBdrFunction(G)
    A.assembleBdr(bc.get("Neumann"), u*g_N.tr() * gismo.expr.nv(G))
  

    solver.compute(A.matrix())

   

    solVector = solver.solve(A.rhs())

    u_sol = A.getSolution(u, solVector)


    l2err = math.sqrt(ev.integral((u_ex - u_sol).sqNorm() * gismo.expr.meas(G)))
    h1err = math.sqrt(ev.integral((gismo.expr.igrad(u_ex) - gismo.expr.igrad(u_sol ,G)).sqNorm() * gismo.expr.meas(G)))

    print(f"L2 error: {l2err}\nH1 error: {h1err}")


