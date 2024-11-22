/** @file composed_domain_L2.cpp

    @brief Tutorial on how to use expression assembler to solve the Poisson equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

//! [Include namespace]
#include <gismo.h>
#include <gsAssembler/gsExprEvaluator.h>
#include <gsCore/gsComposedFunction.h>
#include <gsNurbs/gsSquareDomain.h>

using namespace gismo;
using namespace expr;
//! [Include namespace]

// Integrates one basis function using the expression evaluator
template <short_t DIM, class T>
T testEvaluator(const gsTensorBSpline<DIM,T> & geometry,
                const gsFunction<T>          & composition,
                const gsOptionList           & options,
                const index_t                & idxI,
                const index_t                & idxJ)
{
    gsComposedGeometry<T> cgeom(composition,geometry);
    gsMultiBasis<T> ib(geometry.basis());

    gsMultiPatch<T> mp(geometry);
    gsMultiBasis<T> mb(mp);

    gsExprEvaluator<T> ev;
    ev.setIntegrationElements(ib);
    // Functions i and j
    gsBasisFun<T> Fi = mb.basis(0).function(idxI);
    gsBasisFun<T> Fj = mb.basis(0).function(idxJ);
    ev.options().update(options);
    auto G = ev.getMap(mp);
    auto fi = ev.getVariable(Fi);
    auto fj = ev.getVariable(Fj);
    gsDebugVar(ev.integral(fi*fj*meas(G)));

    gsMultiPatch<T> cmp(cgeom);
    gsMultiBasis<T> cmb(cmp);

    gsExprEvaluator<T> cev;
    cev.setIntegrationElements(ib);
    gsBasisFun<T> cFi = cmb.basis(0).function(idxI);
    gsBasisFun<T> cFj = cmb.basis(0).function(idxJ);
    cev.options().update(options);
    auto cG = cev.getMap(cmp);
    auto cfi = cev.getVariable(cFi);
    auto cfj = cev.getVariable(cFj);
    gsDebugVar(cev.integral(cfi*cfj*meas(cG)));

    return math::pow(ev.integral(fi*fj*meas(G)) - cev.integral(cfi*cfj*meas(cG)),2);
}

// Integrates all basis functions using the expression assembler
template <short_t DIM, class T>
T testIntegrator(const gsTensorBSpline<DIM,T> & geometry,
                 const gsFunction<T>          & composition,
                 const gsOptionList           & options)
{
    gsComposedGeometry<T> cgeom(composition,geometry);
    gsMultiBasis<T> ib(geometry.basis());

    gsMultiPatch<T> mp(geometry);
    gsMultiBasis<T> mb(mp);
    gsExprAssembler<T> A(1,1);
    A.setIntegrationElements(ib);
    A.options().update(options);
    auto G = A.getMap(mp);
    auto u = A.getSpace(mb);
    u.setup();
    A.initSystem();
    A.assemble(u*u.tr()*meas(G));
    gsDebugVar(A.matrix().toDense());

    gsMultiPatch<T> cmp(cgeom);
    gsMultiBasis<T> cmb(cmp);

    gsExprAssembler<T> cA(1,1);
    cA.setIntegrationElements(ib);
    cA.options().update(options);
    auto cG = cA.getMap(cmp);
    auto cu = cA.getSpace(cmb);
    cu.setup();
    cA.initSystem();
    cA.assemble(cu*cu.tr()*meas(cG));
    gsDebugVar(cA.matrix().toDense());

    return (A.matrix()-cA.matrix()).norm();
}

int main(int argc, char *argv[])
{
    index_t degree = 1;
    index_t interior = 0;
    index_t multiplicity = 1;
    index_t nGauss = -1;
    index_t bIndexI = 0;
    index_t bIndexJ = 0;
    std::string X = "x^1";
    std::string Y = "y^1";

    gsCmdLine cmd("");
    cmd.addInt( "p", "degree", "Degree of the basis", degree );
    cmd.addInt( "i", "interior", "Interior knots in the basis", interior );
    cmd.addInt( "m", "multiplicity", "Interior knot multiplicity in the basis", multiplicity );
    cmd.addInt( "G", "numGauss", "Number of Gauss points (default = degree+1)", nGauss );
    cmd.addInt( "I", "bIdxI", "Basis function index of the basis function to test", bIndexI );
    cmd.addInt( "J", "bIdxJ", "Basis function index of the basis function to test", bIndexJ );
    cmd.addString( "X", "xcoord", "X-coordinate of the composition", X );
    cmd.addString( "Y", "ycoord", "Y-coordinate of the composition", Y );

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    nGauss = (nGauss==-1) ? degree + 1 : nGauss;

    gsFunctionExpr<> composition(X,Y,2);
    gsKnotVector<> kv(0,1,interior,degree+1,multiplicity,degree);
    gsDebugVar(kv);
    gsTensorBSplineBasis<2,real_t> tbasis(kv,kv);
    gsMatrix<> coefs = tbasis.anchors().transpose();
    gsTensorBSpline<2,real_t> tspline(tbasis,coefs);

    gsOptionList opt;
    opt.addReal("quA", "Number of quadrature points: quA*deg + quB", 0  );
    opt.addInt ("quB", "Number of quadrature points: quA*deg + quB", nGauss    );
    opt.addInt ("plot.npts", "Number of sampling points for plotting", 3000 );
    opt.addSwitch("plot.elements", "Include the element mesh in plot (when applicable)", false);
    opt.addSwitch("flipSide", "Flip side of interface where evaluation is performed.", false);

    testEvaluator(tspline,composition,opt,bIndexI,bIndexJ);
    testIntegrator(tspline,composition,opt);

    return EXIT_SUCCESS;

}// end main
