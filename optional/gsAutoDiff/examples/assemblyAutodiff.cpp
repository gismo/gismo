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

    // Define a function
    real_t theta = 0.5;
    gsConstantFunction<> fun(theta*theta,1);

    // Define a domain
    gsBSpline<>         geom = *gsNurbsCreator<>::BSplineUnitInterval(1);
    gsBSplineBasis<>    basis= geom.basis();
    basis.degreeElevate(numElevate);
    for (index_t i=0; i!=numHRef; i++)
        basis.uniformRefine();

    gsMatrix<> rhs = computeRhs(geom,basis,fun);
    gsInfo<<rhs.transpose()<<"\n";

    //////////////////////////////////////////
    // Forward-mode AD example
    //////////////////////////////////////////
    dual theta_f = dual(theta);
    autodiff::detail::seed<1>(theta_f, 1.0);

    gsConstantFunction<dual> fun_f(theta_f*theta_f,1);
    gsBSpline<dual>         geom_f = *gsNurbsCreator<dual>::BSplineUnitInterval(1);
    gsBSplineBasis<dual>    basis_f= geom_f.basis();
    basis_f.degreeElevate(numElevate);
    for (index_t i=0; i!=numHRef; i++)
        basis_f.uniformRefine();

    gsMatrix<dual> rhs_f = computeRhs(geom_f,basis_f,fun_f);
    gsInfo<<rhs_f.transpose().cast<double>()<<"\n";
    gsInfo<<"dRHS/dtheta (forward-mode):\n";
    for (index_t i=0; i<rhs_f.rows(); ++i)
    {
        gsInfo<<autodiff::detail::derivative<1>(rhs_f(i,0))<<" ";
    }
    gsInfo<<"\n\n";

    //////////////////////////////////////////
    // Reverse-mode AD example
    //////////////////////////////////////////
    var theta_v = var(theta);
    gsConstantFunction<var> fun_v(theta_v*theta_v,1);
    gsBSpline<var>         geom_v = *gsNurbsCreator<var>::BSplineUnitInterval(1);
    gsBSplineBasis<var>    basis_v= geom_v.basis();
    basis_v.degreeElevate(numElevate);
    for (index_t i=0; i!=numHRef; i++)
        basis_v.uniformRefine();

    gsMatrix<var> rhs_v = computeRhs(geom_v,basis_v,fun_v);
    gsInfo<<rhs_v.transpose().cast<double>()<<"\n";
    gsInfo<<"dRHS/dtheta (reverse-mode):\n";
    for (index_t i=0; i<rhs_v.rows(); ++i)
    {
        auto grads = autodiff::derivatives(rhs_v(i,0), autodiff::reverse::detail::wrt(theta_v));
        gsInfo<<grads[0]<<" ";
    }

    return 0;
}
