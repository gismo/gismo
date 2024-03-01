"""
    @file Heat equation example

    @brief Solve the heat equation using the expression assembler in Python. Needs the gismo cppyy bindings.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Scholz
"""


from gismo_cppyy import gismo
import matplotlib.pyplot as plt
import numpy as np

def plotField(field, evaluator, time, ax, n=10):
    corners = gismo.gsMatrix["real_t"].fromnumpy(np.array([[0., 1.], [0.,1.]]))
    numpoints = gismo.gsVector["int"](2)
    numpoints[0] = n
    numpoints[1] = n
    gismogrid = gismo.gsGridIterator["double", gismo.CUBE, 2](corners, numpoints)

    evaluator.eval(field, gismogrid, 0)
    evalpoints = gismo.gsVector["real_t"](evaluator.allValues()).tonumpy()

    ax.set_title(f"Time: {time}s", loc="left")
    ax.imshow(evalpoints.reshape(n, n), interpolation="bilinear", cmap="coolwarm", vmin=0.0, vmax=.5)


def main():
    numRefine = 5
    theta = 0.5  # 0.0 explicit Euler, 1.0 implicit Euler, 0.5 Crank Nicolson
    endTime = 1.
    numSteps = 400

    plotresolution = 100

    Dt = endTime / numSteps  # Time step

    multipatch = gismo.gsMultiPatch["real_t"](gismo.gsNurbsCreator["real_t"].BSplineSquareDeg(2))

    basis = gismo.gsMultiBasis["real_t"](multipatch, True)

    boundarycondition = gismo.gsBoundaryConditions["real_t"]()
    neumanncond = gismo.gsFunctionExpr["real_t"] ("-1", 2)
    dirichletcond = gismo.gsFunctionExpr["real_t"] ("0", 2)

    boundarycondition.addCondition(0, gismo.boundary.west, gismo.condition_type.neumann, neumanncond)
    boundarycondition.addCondition(0, gismo.boundary.east, gismo.condition_type.dirichlet, dirichletcond)
    boundarycondition.addCondition(0, gismo.boundary.north, gismo.condition_type.dirichlet, dirichletcond)
    boundarycondition.addCondition(0, gismo.boundary.south, gismo.condition_type.dirichlet, dirichletcond)
    boundarycondition.setGeoMap(multipatch)

    print(boundarycondition)

    for _ in range(numRefine):
        basis.uniformRefine()

    # Assemble system matrix
    A_stiff = gismo.gsExprAssembler["real_t"]()
    A_stiff.setIntegrationElements(basis)
    G_stiff = A_stiff.getMap(multipatch)
    u_stiff = A_stiff.getSpace(basis)
    g_N_stiff = A_stiff.getBdrFunction(G_stiff)
    u_stiff.setup(boundarycondition, gismo.dirichlet.l2Projection, 0)
    A_stiff.initSystem()
    A_stiff.assemble(gismo.expr.igrad(u_stiff, G_stiff) * gismo.expr.igrad(u_stiff, G_stiff).tr() * gismo.expr.meas(G_stiff))
    A_stiff.assembleBdr(boundarycondition.get("Neumann"), u_stiff * g_N_stiff.tr() * gismo.expr.nv(G_stiff))

    systemMatrix = A_stiff.matrix()

    # Assemble mass matrix
    A_mass = gismo.gsExprAssembler["real_t"]()
    A_mass.setIntegrationElements(basis)
    G_mass = A_mass.getMap(multipatch)
    u_mass = A_mass.getSpace(basis)
    u_mass.setup(boundarycondition, gismo.dirichlet.l2Projection, 0)
    A_mass.initSystem()
    A_mass.assemble(u_mass*u_mass.tr() * gismo.expr.meas(G_mass))

    massMatrix = A_mass.matrix()

    solver = gismo.gsSparseSolver["real_t"].LU()
    solvector = gismo.gsMatrix["real_t"]()
    solvector.setZero(A_stiff.numDofs(), 1)

    ev = gismo.gsExprEvaluator["real_t"](A_stiff)

    # Initial plot
    plt.ion()
    fig = plt.figure()
    ax = fig.add_subplot(111)

    u_sol = A_stiff.getSolution(u_stiff, solvector)
    plotField(field=u_sol, evaluator=ev, time=0.0, ax=ax, n=plotresolution)
    fig.canvas.draw()
    fig.canvas.flush_events()

    solver.compute(massMatrix + systemMatrix * Dt * theta)
    for step in range(numSteps):
        solvector = solver.solve(A_stiff.rhs() * Dt + (massMatrix - systemMatrix * Dt*(1-theta)) * solvector)

        u_sol = A_stiff.getSolution(u_stiff, solvector)
        plotField(field=u_sol, evaluator=ev, time=step*Dt, ax=ax, n=plotresolution)
        fig.canvas.draw()
        fig.canvas.flush_events()
        plt.pause(1/24)
    plt.show()


if __name__ == "__main__":
    main()