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

template<typename T> using Real       = codi::RealReverseGen<T>;
template<typename T> using Tape       = typename Real<T>::Tape;
template<typename T> using Identifier = typename Real<T>::Identifier;
template<typename T> using RealBase   = typename Real<T>::Real;

// Implementation of the system solver for CoDiPack
template<class T>
void solveSystem_rev(Tape<T>                                                            *tape,
                     void                                                               *d,
                     codi::VectorAccessInterface<typename Tape<T>::Real, Identifier<T>> *va)
{
    codi::ExternalFunctionUserData* data = (codi::ExternalFunctionUserData*)d;
    gsSparseMatrix<Real<T>>* matrix;
    gsMatrix<Real<T>>* rhs;
    gsMatrix<Real<T>>* sol;

    // Step 4: Get data in the same order as it was added.
    data->getData(matrix);
    data->getData(rhs);
    data->getData(sol);

    typename gsSparseSolver<Real<T>>::LU solver;

    gsMatrix<Real<T>> solAdj(*sol);
    gsMatrix<Real<T>> rhsAdj(*sol);
    for(index_t i = 0; i < sol->size(); ++i) {
        solAdj[i] = (*sol)[i].getGradient();
    }

    solver.compute(matrix->transpose());
    rhsAdj = solver.solve(solAdj);

    for(index_t i = 0; i < sol->size(); ++i) {
        auto index = (*rhs)[i].getIdentifier();
        tape->gradient(index) += rhsAdj[i].getValue();
    }
    for (int e=0; e<matrix->outerSize(); ++e) {
        for (typename gsSparseMatrix<Real<T>>::InnerIterator it(*matrix,e); it; ++it) {
            index_t k = it.row();
            index_t l = it.col();
            Real<T>& temp1 = matrix->at(k,l);
            tape->gradient(temp1.getIdentifier()) += -rhsAdj[l].getValue() * (*sol)[k].getValue();
        }
    }
}

// Implementation of the system deleter for CoDiPack
template<class T>
void solveSystem_del(Tape<T>* tape, void* d) {
  codi::ExternalFunctionUserData* data = (codi::ExternalFunctionUserData*)d;
  
  // Step 5: Delete the data
  delete data;
}

// Implementation of the system solver for CoDiPack
template<class T>
void solveSystem(typename gsSparseSolver<Real<T>>::LU &solver,
                 const gsSparseMatrix<Real<T>>        &matrix,
                 const gsMatrix<Real<T>>              &rhs,
                 gsMatrix<Real<T>>                    &sol) {

    Tape<T>& tape = Real<T>::getTape();
    tape.setPassive();

    // Step 1: Create the data object
    codi::ExternalFunctionUserData* data = new codi::ExternalFunctionUserData();

    // Step 2: Add data
    data->addData(&matrix);
    data->addData(&rhs);
    
    solver.compute( matrix );
    sol = solver.solve( rhs );
    data->addData(&sol);

    // Step 3: Add the external function with the data
    tape.pushExternalFunction(codi::ExternalFunction<Tape<T>>::create(solveSystem_rev<T>, data, solveSystem_del<T>));
    tape.setActive();
    for(index_t i = 0; i < sol.size(); ++i) {
        tape.registerInput(sol[i]);
    }
}

int main(int argc, char* argv[])
{
    // Input options
    int numElevate  = 0;
    int numRefine   = 0;
    int basisDegree = 0;
    bool effSolve   = false;

    gsCmdLine cmd("Testing the CoDiPack extension.");
    cmd.addInt("r","refine",
               "Number of dyadic h-refinement (bisection) steps to perform before solving",
               numRefine);
    cmd.addInt("p","degree",
               "Degree of the basis functions to use for solving (will elevate or reduce the input)",
               basisDegree);
    cmd.addInt("e","elevate",
               "Number of degree elevation steps to perform on the geometry's basis before solving",
               numElevate);
    cmd.addSwitch("effSolve", "Solve the system efficiently", effSolve);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    Real<real_t>::Tape& tape = Real<real_t>::getTape();

    Real<real_t> a = 3.0;
    Real<real_t> b = 2.0;

    tape.setActive();
    tape.registerInput(a);
    tape.registerInput(b);

    gsMultiPatch<Real<real_t>> mp( *gsNurbsCreator<Real<real_t>>::BSplineRectangle(0.0,0.0,a,b) );

    gsFunctionExpr<Real<real_t>>  f("0.0", 2);
    gsFunctionExpr<Real<real_t>>  g("1.0", 2);
    gsFunctionExpr<Real<real_t>>  coeff_A("1.0","0","0","1.0", 2);
    gsFunctionExpr<Real<real_t>>  coeff_b("0.2","0.4", 2);
    gsFunctionExpr<Real<real_t>>  coeff_c("0", 2);

    gsBoundaryConditions<Real<real_t>> BCs;
    for (gsMultiPatch<Real<real_t>>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
    {
        BCs.addCondition( *bit, condition_type::dirichlet, &g);
    }

    gsMultiBasis<Real<real_t>> bases(mp);
    if (basisDegree)
        bases.setDegree(basisDegree);
    else if (numElevate)
        bases.degreeElevate(numElevate);
    for (int i=0; i<numRefine ; i++)
        bases.uniformRefine();

    gsSparseSolver<Real<real_t>>::LU solver;
    gsCDRAssembler<Real<real_t>> galerkin(mp, bases, BCs, f, coeff_A, coeff_b, coeff_c,
                                          dirichlet::elimination, iFace::glue, stabilizerCDR::none);
    galerkin.assemble();
    gsMatrix<Real<real_t>> solVector;

    std::cout << "\n\nTape statistics before solving the system:\n\n";
    tape.printStatistics();

    if (effSolve)
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

    gsField<Real<real_t>> sol = galerkin.constructSolution(solVector);

    gsConstantFunction<Real<real_t>> zero (0.0,2);
    Real<real_t> result = sol.distanceL2(zero,true);

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
