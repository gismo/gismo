/** @file codipack_example

    @brief Tutorial on algorithmic differentiation in reverse mode
    with the use of CoDiPack with memory-efficient solution of adjoint
    linear system

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Jaeschke
*/

#include <iostream>

#include <gismo.h>

using namespace gismo;

// Implementation of the system solver for CoDiPack
void solveSystem_b(codi::RealReverse::TapeType         *tape,
                   codi::DataStore                     *data,
                   codi::AdjointInterface<double, int> *interface) {

    gsSparseMatrix<codi::RealReverse>* matrix;
    gsMatrix<codi::RealReverse>* rhs;
    gsMatrix<codi::RealReverse>* sol;

    data->getData(matrix);
    data->getData(rhs);
    data->getData(sol);

    gsSparseSolver<codi::RealReverse>::LU solver;

    gsMatrix<codi::RealReverse> solAdj(*sol);
    gsMatrix<codi::RealReverse> rhsAdj(*sol);
    for(index_t i = 0; i < sol->size(); ++i) {
        solAdj[i] = (*sol)[i].getGradient();
    }

    solver.compute(matrix->transpose());
    rhsAdj = solver.solve(solAdj);

    for(index_t i = 0; i < sol->size(); ++i) {
        auto index = (*rhs)[i].getGradientData();
        tape->gradient(index) += rhsAdj[i].getValue();
    }
    for (int e=0; e<matrix->outerSize(); ++e) {
        for (gsSparseMatrix<codi::RealReverse>::InnerIterator it(*matrix,e); it; ++it) {
            int k = it.row();
            int l = it.col();
            codi::RealReverse& temp1 = matrix->at(k,l);
            tape->gradient(temp1.getGradientData()) += -rhsAdj[l].getValue() * (*sol)[k].getValue();
        }
    }
}

// Implementation of the system deleter for CoDiPack
void solveSystem_delete(codi::RealReverse::TapeType* tape, codi::DataStore* data) {}

// Implementation of the system solver for adDSL
void solveSystem(gsSparseSolver<codi::RealReverse>::LU   &solver,
                 const gsSparseMatrix<codi::RealReverse> &matrix,
                 const gsMatrix<codi::RealReverse>       &rhs,
                 gsMatrix<codi::RealReverse>             &sol) {

    codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
    tape.setPassive();

    codi::DataStore* dataHandler = new codi::DataStore();
    dataHandler->addData(&matrix);
    dataHandler->addData(&rhs);
    solver.compute( matrix );
    sol = solver.solve( rhs );

    dataHandler->addData(&sol);

    tape.pushExternalFunction(&solveSystem_b, dataHandler, &solveSystem_delete);
    tape.setActive();
    for(index_t i = 0; i < sol.size(); ++i) {
        tape.registerInput(sol[i]);
    }
}

int main(int argc, char* argv[])
{
    // Input options
    int numElevate  = 0;
    int numHref     = 0;
    int basisDegree = 0;
    bool EffSolv    = false;

    gsCmdLine cmd("Testing compressible Euler problem.");
    cmd.addInt("r","hRefine",
               "Number of dyadic h-refinement (bisection) steps to perform before solving",
               numHref);
    cmd.addInt("p","degree",
               "Degree of the basis functions to use for solving (will elevate or reduce the input)",
               basisDegree);
    cmd.addInt("e","degreeElevation",
               "Number of degree elevation steps to perform on the Geometry's basis before solving",
               numElevate);
    cmd.addSwitch("effSolv", "Solve the system efficiently", EffSolv);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();

    codi::RealReverse a = 3.0;
    codi::RealReverse b = 2.0;

    tape.setActive();
    tape.registerInput(a);
    tape.registerInput(b);

    gsMultiPatch<codi::RealReverse> mp( *gsNurbsCreator<codi::RealReverse>::BSplineRectangle(0.0,0.0,a,b) );

    gsFunctionExpr<codi::RealReverse>  f("0.0", 2);
    gsFunctionExpr<codi::RealReverse>  g("1.0", 2);
    gsFunctionExpr<codi::RealReverse>  coeff_A("1.0","0","0","1.0", 2);
    gsFunctionExpr<codi::RealReverse>  coeff_b("0.2","0.4", 2);
    gsFunctionExpr<codi::RealReverse>  coeff_c("0", 2);

    gsBoundaryConditions<codi::RealReverse> BCs;
    for (gsMultiPatch<codi::RealReverse>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
    {
        BCs.addCondition( *bit, condition_type::dirichlet, &g);
    }

    gsMultiBasis<codi::RealReverse> bases(mp);
    if (basisDegree)
        bases.setDegree(basisDegree);
    else if (numElevate)
        bases.degreeElevate(numElevate);
    for (int i=0; i<numHref ; i++)
        bases.uniformRefine();

    gsSparseSolver<codi::RealReverse>::LU solver;
    gsCDRAssembler<codi::RealReverse> galerkin(mp, bases, BCs, f, coeff_A, coeff_b, coeff_c,
                                               dirichlet::elimination, iFace::glue, stabilizerCDR::none);
    galerkin.assemble();
    gsMatrix<codi::RealReverse> solVector;

    std::cout << "\n\nTape statistics before solving the system:\n\n";
    tape.printStatistics();

    if (EffSolv)
    {
        // Efficient way of solving the system
        solveSystem(solver, galerkin.matrix(), galerkin.rhs(), solVector);
    }
    else
    {
        // Inefficient way of solving the system
        solver.compute(galerkin.matrix());
        solVector = solver.solve(galerkin.rhs());
    }

    std::cout << "\n\nTape statistics after solving the system:\n\n";
    tape.printStatistics();

    gsField<codi::RealReverse> sol = galerkin.constructSolution(solVector);

    gsConstantFunction<codi::RealReverse> zero (0.0,2);
    codi::RealReverse result = sol.distanceL2(zero,true);

    result = result * result;
    tape.registerOutput(result);
    tape.setPassive();

    std::cout << "\n\nTape statistics at the end:\n\n";
    tape.printStatistics();

    result.setGradient(1.0);

    tape.evaluate();
    std::cout << "\n\nThis code computes the derivatives of area of a 3x2 rectangle with respect to lengths of its sides:\n\n";

    std::cout << a.getGradient() << std::endl
              << b.getGradient() << std::endl;

    return 0;
}

/*
int main(int argc, char* argv[])
{
    gsInfo << "Compile with -DGISMO_WITH_CODIPACK=ON -DGISMO_BUILD_LIB=OFF\n";
    return 0;
}
*/
