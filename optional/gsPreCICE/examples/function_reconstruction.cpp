/** @file flow-over-heated-plate.cpp

    @brief Heat equation participant for the PreCICE example "flow over heated plate"

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 0;
    index_t numElevate = 0;

    gsCmdLine cmd("Flow over heated plate for PreCICE.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
  
    gsFunctionExpr<> fun("sin(x*pi)*cos(y*pi)*exp(x*y)",2);
    gsMultiPatch<> mp;
    mp.addPatch(gsNurbsCreator<>::BSplineSquare());

    mp.degreeElevate(numElevate);
    for (index_t k=0; k!=numRefine; k++)
        mp.uniformRefine();


    gsKnotVector<>::knotContainer kc0 = dynamic_cast<gsBSplineBasis<real_t> *>(&mp.basis(0).component(0))->knots().unique();
    gsKnotVector<>::knotContainer kc1 = dynamic_cast<gsBSplineBasis<real_t> *>(&mp.basis(0).component(1))->knots().unique();
    gsKnotVector<> kv0(kc0,0);
    gsKnotVector<> kv1(kc1,0);

    // gsTensorBSplineBasis<2,real_t> tbasis(kv0,kv1);
    const gsBasis<> & tbasis = mp.basis(0);

    gsMultiBasis<> mb(mp);

    gsField<> field(mp,fun);
    gsWriteParaview(field,"field",5000);

    gsExprEvaluator<> ev;
    ev.setIntegrationElements(mb);
    auto G = ev.getMap(mp);
    auto f = ev.getVariable(fun,G);
    typedef typename gsExprHelper<real_t>::element     element;     ///< Current element
    element el= ev.getElement();

    gsDebugVar(ev.integral(f*meas(G)));
    ev.integralElWise(f*meas(G));
    gsDebugVar(gsAsConstVector<real_t>(ev.elementwise()));
    gsDebugVar(gsAsConstVector<real_t>(ev.elementwise()).sum());

    gsAsConstVector<real_t> F(ev.elementwise());

    gsMultiPatch<> mp0;
    mp0.addPatch(tbasis.makeGeometry(gsAsConstVector<real_t>(ev.elementwise())));
    gsWriteParaview(mp0,"geom",5000);

    // gsMatrix<> points = mp.basis(0).anchors();
    // gsSparseMatrix<> C = mp.basis(0).collocationMatrix(points);
    // gsMatrix<> F = mp0.piece(0).eval(points).transpose();
    // gsMultiPatch<> mpnew;

    // gsDebugVar(mp.basis(0));
    // mpnew.addPatch(mp.basis(0).makeGeometry(C.transpose()*F));
    // gsWriteParaview(mpnew,"geomnew",5000);

    // auto fsmooth = ev.getVariable(mpnew,G);
    // ev.integralElWise(fsmooth*meas(G) / el.area());
    // gsDebugVar((C.transpose() * F));
    // gsDebugVar(gsAsConstVector<real_t>(ev.elementwise()));

    // gsDebugVar((C.transpose() * F).sum());

    // ev.writeParaview(fsmooth/el.area(),G,"pressure");



    ///////////////////////////////////////////////////////////

    // gsOptionList opt;
    // opt.addReal("quA", "Number of quadrature points: quA*deg + quB; For patchRule: Regularity of the target space", 1.0  );
    // opt.addInt ("quB", "Number of quadrature points: quA*deg + quB; For patchRule: Degree of the target space", 1    );
    // opt.addInt ("quRule", "Quadrature rule used (1) Gauss-Legendre; (2) Gauss-Lobatto; (3) Patch-Rule",1);
    // opt.addSwitch("overInt", "Apply over-integration on boundary elements or not?", false);


    // gsSparseMatrix<> result( mp.basis(0).numElements(), mp.basis(0).size() );
    gsBasis<>::domainIter domIt = tbasis.makeDomainIterator();

    // typename gsQuadRule<>::uPtr QuRule; // Quadrature rule
    // QuRule = gsQuadrature::getPtr(mp.basis(0), opt);
    // gsVector<> quWeights; // quadrature weights
    // gsMatrix<> quPoints; // quadrature weights

    // gsMatrix<index_t> actives;
    // gsMatrix<> vals;
    // gsVector<> integrals;

    gsMatrix<> centerPoints(2,mp.basis(0).numElements());


    // gsMatrix<> evals(mp.basis(0).numElements(),1);

    // gsDebugVar(mp.basis(0).numElements());

    // gsDebugVar(mp.basis(0).size());

// #pragma omp parallel
// {
// #   ifdef _OPENMP
//     const int tid = omp_get_thread_num();
//     const int nt  = omp_get_num_threads();
// #   endif

    // Start iteration over elements of patchInd
// #   ifdef _OPENMP
//     for ( domIt->next(tid); domIt->good(); domIt->next(nt) )
// #   else
     for (; domIt->good(); domIt->next() )
// #   endif
     {

        // Map the Quadrature rule to the element
        // QuRule->mapTo( domIt->lowerCorner(), domIt->upperCorner(),
        //                quPoints, quWeights);

        // mp.basis(0).eval_into(quPoints,vals);
        // mp.basis(0).active_into(domIt->centerPoint(),actives);

        centerPoints.col(domIt->id()) = domIt->centerPoint();

//         integrals = vals * quWeights;

//         for (index_t k=0; k!=actives.rows(); k++)
// #           pragma omp critical (int_colloc)
//             {
//                 // Element*basis function
//                 result(domIt->id(),actives(k,0)) = integrals.at(k);
//             }
     }
//}

    gsMatrix<> RHS = mp0.piece(0).eval(centerPoints).transpose();


    gsSparseMatrix<> result = mp.basis(0).collocationMatrixIntegrated();
    // gsDebugVar(RHS.transpose());
    /*


    /// Next steps: Solve Cx = F. If not square, use pseudo inverse
    // typename gsMatrix<T>::BCSVD bcSvd = derivs.jacobiSvd(gsEigen::ComputeFullV);

    // gsSparseSolver<>::CGDiagonal solver;
    // gsSparseMatrix<> A = result.transpose()*result;

    // gsDebugVar(A.rows());
    // gsDebugVar(A.cols());

    // gsMatrix<> B = result.transpose()*RHS;
    // solver.compute(A);
    // gsMatrix<> X = solver.solve(B);

    // gsDebugVar(X);

    // 

    gsMultiPatch<> mpnew;

    // TODO: Rewrite the collocation solver for the SVD problem
    
    gsMatrix<> dense = result.toDense();
    gsMatrix<>::JacobiSVD jacobiSvd(dense,gsEigen::ComputeFullU | gsEigen::ComputeFullV);
    gsMatrix<> X = jacobiSvd.solve(RHS);


    mpnew.addPatch(mp.basis(0).makeGeometry(X));
    gsWriteParaview(mpnew,"geomnew",5000);

    auto FF = ev.getVariable(mpnew,G);

    ev.integralElWise(FF*meas(G));
    gsDebugVar(gsAsConstVector<real_t>(ev.elementwise()));
    gsDebugVar(gsAsConstVector<real_t>(ev.elementwise()).sum());




    // A.template bdcSvd<Eigen::ComputeThinU | Eigen::ComputeThinV>().solve(b)

    */

    return  EXIT_SUCCESS;
}
