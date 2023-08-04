#!/usr/bin/python

""""
    @file Poisson example

    @brief Solve the Poisson equation in Python

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Scholz
"""

import os, sys
gismo_path=os.path.join(os.path.dirname(__file__), "../cmake-build-debug/lib")
sys.path.append(gismo_path)

import pygismo as gs
import scipy.sparse
import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

numRefine = 3
degree = 3
plot = True

# Right hand side
exactFunc = gs.core.gsFunctionExpr("cos(2*pi*x) + cos(2*pi*y) + 2", 2)
rhs = gs.core.gsFunctionExpr("4*pi^2*cos(2*pi*x) + 4*pi^2*cos(2*pi*y)", 2)

# Make domain
#domain = gs.nurbs.gsNurbsCreator.BSplineSquare(1,0,0)
domain = gs.nurbs.gsNurbsCreator.BSplineQuarterAnnulus(2)

# Multipatch
mp = gs.core.gsMultiPatch(domain)

# Get basis
basis = gs.core.gsMultiBasis(mp)
# Refine and elevate degree
basis.uniformRefine(numRefine, 1, -1)
basis.degreeElevate(max(0,degree - domain.degree(0)), -1)

#Make boundary conditions
bc = gs.pde.gsBoundaryConditions()
bc.addCondition(0, gs.core.side.west, gs.pde.bctype.dirichlet, exactFunc, 0, False, -1)
bc.addCondition(0, gs.core.side.east, gs.pde.bctype.dirichlet, exactFunc, 0, False, -1)
bc.addCondition(0, gs.core.side.south, gs.pde.bctype.dirichlet, exactFunc, 0, False, -1)
bc.addCondition(0, gs.core.side.north, gs.pde.bctype.dirichlet, exactFunc, 0, False, -1)

# Make assembler
assembler = gs.assembler.gsPoissonAssembler(mp,
                                            basis,
                                            bc,
                                            rhs,
                                            gs.assembler.dirichletStrategy.elimination,
                                            gs.assembler.iFaceStrategy.conforming)

# Assemble
assembler.assemble()

# Get matrix and rhs
matrix = assembler.matrix()
rhs = assembler.rhs()

# Solve
sol = scipy.sparse.linalg.spsolve(matrix, rhs)

# Create solution and compute errors
solutionField = assembler.constructSolution(sol, 0)
L2Error = solutionField.distanceL2(exactFunc, False, 1000)
H1SeminormError = solutionField.distanceH1(exactFunc, False, 1000)
print(f"L2 error: {L2Error}\nH1 error: {L2Error + H1SeminormError}")


if plot:
    # Plot the solution
    numsamples = 30
    x = np.linspace(0, 1, numsamples, True)
    y = np.linspace(0, 1, numsamples, True)
    X, Y = np.meshgrid(x,y)

    evalParams = np.concatenate([X.reshape(1,-1), Y.reshape(1,-1)], axis=0)

    evalPoints = solutionField.value(evalParams, 0)
    evalDomain = solutionField.point(evalParams, 0)


    plotX = evalDomain[0,:].reshape(numsamples, numsamples)
    plotY = evalDomain[1,:].reshape(numsamples,numsamples)
    plotZ = evalPoints.reshape(numsamples,numsamples)

    # Evaluation of the exact function
    # exactPoints = exactFunc.eval(evalDomain)
    # plotZexact = exactPoints.reshape(numsamples, numsamples)

    ax = plt.subplot(111, projection='3d')
    ax.plot_surface(plotX, plotY, plotZ)
    ax.plot_surface(plotX, plotY, np.zeros_like(plotX))
    plt.show()
