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
#include <gsCoDiPack/gsCoDiPack.h>

using namespace gismo;

// Implementation of the system solver for CoDiPack
template<class T>
void solveSystem_b(codi::RealReverseGen<T>::TapeType *tape,
                   codi::DataStore                     *data,
                   codi::AdjointInterface<double, int> *interface) {

    gsSparseMatrix<codi::RealReverseGen<T>>* matrix;
    gsMatrix<codi::RealReverseGen<T>>* rhs;
    gsMatrix<codi::RealReverseGen<T>>* sol;

    data->getData(matrix);
    data->getData(rhs);
    data->getData(sol);

    gsSparseSolver<codi::RealReverseGen<T>>::LU solver;

    gsMatrix<codi::RealReverseGen<T>> solAdj(*sol);
    gsMatrix<codi::RealReverseGen<T>> rhsAdj(*sol);
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
        for (gsSparseMatrix<codi::RealReverseGen<T>>::InnerIterator it(*matrix,e); it; ++it) {
            int k = it.row();
            int l = it.col();
            codi::RealReverseGen<T>& temp1 = matrix->at(k,l);
            tape->gradient(temp1.getGradientData()) += -rhsAdj[l].getValue() * (*sol)[k].getValue();
        }
    }
}

// Implementation of the system deleter for CoDiPack
template<class T>
void solveSystem_delete(codi::RealReverseGen<T>::TapeType* tape, codi::DataStore* data) {}

// Implementation of the system solver for CoDiPack
template<class T>
void solveSystem(gsSparseSolver<codi::RealReverseGen<T>>::LU   &solver,
                 const gsSparseMatrix<codi::RealReverseGen<T>> &matrix,
                 const gsMatrix<codi::RealReverseGen<T>>       &rhs,
                 gsMatrix<codi::RealReverseGen<T>>             &sol) {

    codi::RealReverseGen<T>::TapeType& tape = codi::RealReverseGen<T>::getGlobalTape();
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

    codi::RealReverseGen<real_t>::TapeType& tape = codi::RealReverseGen<real_t>::getGlobalTape();

    codi::RealReverseGen<real_t> a = 3.0;
    codi::RealReverseGen<real_t> b = 2.0;

    tape.setActive();
    tape.registerInput(a);
    tape.registerInput(b);

    gsMultiPatch<codi::RealReverseGen<real_t>> mp( *gsNurbsCreator<codi::RealReverseGen<real_t>>::BSplineRectangle(0.0,0.0,a,b) );

    gsFunctionExpr<codi::RealReverseGen<real_t>>  f("0.0", 2);
    gsFunctionExpr<codi::RealReverseGen<real_t>>  g("1.0", 2);
    gsFunctionExpr<codi::RealReverseGen<real_t>>  coeff_A("1.0","0","0","1.0", 2);
    gsFunctionExpr<codi::RealReverseGen<real_t>>  coeff_b("0.2","0.4", 2);
    gsFunctionExpr<codi::RealReverseGen<real_t>>  coeff_c("0", 2);

    gsBoundaryConditions<codi::RealReverseGen<real_t>> BCs;
    for (gsMultiPatch<codi::RealReverseGen<real_t>>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
    {
        BCs.addCondition( *bit, condition_type::dirichlet, &g);
    }

    gsMultiBasis<codi::RealReverseGen<real_t>> bases(mp);
    if (basisDegree)
        bases.setDegree(basisDegree);
    else if (numElevate)
        bases.degreeElevate(numElevate);
    for (int i=0; i<numHref ; i++)
        bases.uniformRefine();

    gsSparseSolver<codi::RealReverseGen<real_t>>::LU solver;
    gsCDRAssembler<codi::RealReverseGen<real_t>> galerkin(mp, bases, BCs, f, coeff_A, coeff_b, coeff_c,
                                               dirichlet::elimination, iFace::glue, stabilizerCDR::none);
    galerkin.assemble();
    gsMatrix<codi::RealReverseGen<real_t>> solVector;

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

    gsField<codi::RealReverseGen<real_t>> sol = galerkin.constructSolution(solVector);

    gsConstantFunction<codi::RealReverseGen<real_t>> zero (0.0,2);
    codi::RealReverseGen<real_t> result = sol.distanceL2(zero,true);

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
