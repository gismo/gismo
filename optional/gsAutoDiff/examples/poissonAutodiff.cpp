/** @file linearSolveAutoDiff.cpp

    @brief Example demonstrating automatic differentiation through linear solvers

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst

*/

//! [Include namespace]
#include <gismo.h>
#include <gsAutoDiff/gsAutoDiff2.h>
#include <gsAutoDiff/gsAutoDiffUtils.h>
// Include reverse-mode AD with full GISMO matrix support
#include <gsAutoDiff/gsAutoDiffVar.h>
#include <gsAssembler/gsPoissonAssembler.h>

using namespace gismo;
using namespace autodiff;

using namespace gismo;
//! [Include namespace]

using namespace autodiff;

template <class T>
gsMatrix<T> computeRhs(const gsGeometry<T>  & geom,
                       const gsBasis<T>     & basis,
                       const gsFunction<T>  & fun,
                       const gsBoundaryConditions<T>& bc = gsBoundaryConditions<T>() )
{
    gsMultiPatch<T> mp(geom);
    gsMultiBasis<T> mb(basis);
    gsPoissonAssembler<T> assembler(mp, mb, bc,fun);
    assembler.options().setInt("DirichletValues",dirichlet::l2Projection);
    assembler.assemble();
    return assembler.rhs();
}

int main()
{   
    index_t numHRef = 2;
    index_t numElevate = 1;
    // gsCmdLine cmd("AAA");
    // cmd.addReal()

    //////////////////////////////////////////
    // Forward-mode AD example
    //////////////////////////////////////////
    dual theta = dual(0.5);
    autodiff::detail::seed<1>(theta, 1.0);

    gsConstantFunction<dual>fun(theta*theta,1);
    gsBSpline<dual>         geom = *gsNurbsCreator<dual>::BSplineUnitInterval(1);
    gsBSplineBasis<dual>    basis= geom.basis();
    basis.degreeElevate(numElevate);
    for (index_t i=0; i!=numHRef; i++)
        basis.uniformRefine();

    gsMultiPatch<dual> mp(geom);
    gsMultiBasis<dual> mb(basis);

    gsBoundaryConditions<dual> bc;
    bc.addCondition( boundary::west, condition_type::dirichlet, 0 );
    bc.addCondition( boundary::east, condition_type::dirichlet, 0 );
    bc.setGeoMap(mp);

    gsPoissonAssembler<dual> assembler(mp, mb, bc,fun);
    assembler.options().setInt("DirichletValues",dirichlet::l2Projection);
    assembler.assemble();

    gsMatrix<dual> rhs = assembler.rhs();
    gsSparseMatrix<dual> mat = assembler.matrix();

    // Solve with adjoints
    

    gsDebugVar(mat);
    gsDebugVar(rhs);

    gsSparseSolver<dual>::LU solver;
    solver.compute(mat);
    gsMatrix<dual> sol = solver.solve(rhs);

    dual norm = sol.norm();
    
    // COMPUTE DERIVATIVES
    double dnorm_dtheta = autodiff::detail::derivative<1>(norm);
    gsInfo<<"Forward-mode AD result:\n";
    gsInfo<<"Norm of solution: "<<autodiff::detail::derivative<0>(norm)<<"\n";
    gsInfo<<"dNorm/dtheta: "<<dnorm_dtheta<<"\n";


    // gsInfo<<rhs.transpose().cast<double>()<<"\n";
    // gsInfo<<"dRHS/dtheta (forward-mode):\n";
    // for (index_t i=0; i<rhs.rows(); ++i)
    // {
    //     gsInfo<<autodiff::detail::derivative<1>(rhs(i,0))<<" ";
    // }
    // gsInfo<<"\n\n";

    return 0;
}
