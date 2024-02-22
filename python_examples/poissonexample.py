""""
    @file Poisson equation example

    @brief Solve the Poisson equation using the expression assembler in Python. Needs the gismo cppyy bindings.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Scholz
"""

import math
import argparse
import time
from gismo_cppyy import gismo


def main():
    parser = argparse.ArgumentParser(
        prog="poissonexample",
        description="Tutorial on solving a Poisson problem."
    )
    parser.add_argument("-e", "--degreeElevation",
                        default=0,
                        help="Number of degree elevation steps to perform before solving "
                             "(0: equalize degree in all directions)",
                        type=int,
                        dest="numElevate")
    parser.add_argument("-r", "--uniformRefine",
                        default=5,
                        help="Number of Uniform h-refinement loops",
                        type=int,
                        dest="numRefine")
    parser.add_argument("-f", "--file",
                        default="pde/poisson2d_bvp.xml",
                        help="Input XML file",
                        type=str,
                        dest="fn")
    parser.add_argument("--last",
                        dest="last",
                        action="store_const",
                        help="Solve solely for the last level of h-refinement",
                        const=True, default=False)
    parser.add_argument("--plot",
                        dest="plot",
                        action="store_const",
                        help="Create a ParaView visualization file with the solution",
                        const=True, default=False)

    args = parser.parse_args()

    numElevate = args.numElevate
    numRefine = args.numRefine
    fn = args.fn
    last = args.last
    plot = args.plot

    # Read the data from file
    fd = gismo.gsFileData["double"](fn)

    # Get the multipatch geometry
    geometry = gismo.gsMultiPatch['double']()
    fd.getId(0, geometry)

    # Get the right-hand side
    righthandside = gismo.gsFunctionExpr['double']()
    fd.getId(1, righthandside)
    print(f"Source function {righthandside}")

    # Load the boundary conditions
    boundarycondition = gismo.gsBoundaryConditions["double"]()
    fd.getId(2, boundarycondition)
    boundarycondition.setGeoMap(geometry)
    print(f"Boundary conditions:\n{boundarycondition}")

    # Get the exact solution
    exactsolution = gismo.gsFunctionExpr["double"]()
    fd.getId(3, exactsolution)

    # Get options
    Aopt = gismo.gsOptionList()
    fd.getId(4, Aopt)

    # Create basis for solution
    dbasis = gismo.gsMultiBasis["double"](geometry, True)
    dbasis.setDegree(dbasis.maxCwiseDegree() + numElevate)

    # h-refine each basis
    if last:
        for _ in range(numRefine):
            dbasis.uniformRefine()
        numRefine = 0

    print(f"Number of patches: {geometry.nPatches()}, degree for solution: {dbasis.minCwiseDegree()}")

    # Get expression assembler
    A = gismo.gsExprAssembler["double"]()
    A.setOptions(Aopt)
    print(f"Active Options:\n{A.options()}")

    # Elements used for numerical integration
    A.setIntegrationElements(dbasis)

    # Expression evaluator
    ev = gismo.gsExprEvaluator["double"](A)

    # Expressions
    G = A.getMap(geometry)
    u = A.getSpace(dbasis)
    ff = A.getCoeff(righthandside, G)
    u_ex = ev.getVariable(exactsolution, G)

    # Solver
    solver = gismo.gsSparseSolver["double"].CGDiagonal()

    for r in range(numRefine + 1):
        print(f"Refinement step {r}")

        startingtime = time.time()

        if r > 0:
            dbasis.uniformRefine()

        # Set up the space and initialize
        u.setup(boundarycondition, gismo.dirichlet.l2Projection, 0)
        A.initSystem()


        print(f"Number of degrees of freedom: {A.numDofs()}")

        # Assemble the bilinear form
        A.assemble(gismo.expr.igrad(u, G) * gismo.expr.igrad(u, G).tr() * gismo.expr.meas(G),
                   u * ff * gismo.expr.meas(G))
        g_N = A.getBdrFunction(G)
        A.assembleBdr(boundarycondition.get("Neumann"), u * g_N.tr() * gismo.expr.nv(G))

        # Solve
        solver.compute(A.matrix())
        solVector = solver.solve(A.rhs())
        u_sol = A.getSolution(u, solVector)

        # Compute the L2 and H1-errors
        l2error = math.sqrt(ev.integral((u_ex - u_sol).sqNorm() * gismo.expr.meas(G)))
        h1error = math.sqrt(
            ev.integral((gismo.expr.igrad(u_ex) - gismo.expr.igrad(u_sol, G)).sqNorm() * gismo.expr.meas(G)))

        totaltime = time.time() - startingtime

        print(f"L2 error: {l2error}\nH1 error: {h1error}\n"
              f"Total computation time: {totaltime}", flush=True)

    if plot:
        print("Plotting in Paraview...\n")
        collection = gismo.gsParaviewCollection("ParaviewOutput/solution", ev)
        collection.options().setSwitch("plotElements", True)
        collection.options().setInt("plotElements.resolution", 16)
        collection.newTimeStep(geometry)
        collection.addField(u_sol, "numerical solution")
        collection.addField(u_ex, "exact solution")
        collection.saveTimeStep()
        collection.save()
    else:
        print("Done. No output created, re-run with --plot to get a ParaView file containing the solution.")


if __name__ == "__main__":
    main()
