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

    gsTensorBSplineBasis<2,real_t> tbasis(kv0,kv1);


    gsMultiBasis<> mb(mp);

    gsField<> field(mp,fun);
    gsWriteParaview(field,"field",5000);

    gsExprEvaluator<> ev;
    ev.setIntegrationElements(mb);
    auto G = ev.getMap(mp);
    auto f = ev.getVariable(fun,G);
    gsDebugVar(ev.integral(f*meas(G)));
    ev.integralElWise(f*meas(G));
    gsDebugVar(gsAsConstVector<real_t>(ev.elementwise()));
    gsDebugVar(gsAsConstVector<real_t>(ev.elementwise()).sum());

    gsMultiPatch<> mp0;
    mp0.addPatch(tbasis.makeGeometry(gsAsConstVector<real_t>(ev.elementwise())));
    gsWriteParaview(mp0,"geom",5000);

    gsMatrix<> points = mp.basis(0).anchors();
    gsSparseMatrix<> C = mp.basis(0).collocationMatrix(points);
    gsMatrix<> F = mp0.piece(0).eval(points).transpose();
    gsMultiPatch<> mpnew;

    gsDebugVar(mp.basis(0));
    mpnew.addPatch(mp.basis(0).makeGeometry(C.transpose()*F));
    gsWriteParaview(mpnew,"geomnew",5000);

    auto fsmooth = ev.getVariable(mpnew,G);
    ev.integralElWise(fsmooth*meas(G));
    gsDebugVar((C.transpose() * F));
    gsDebugVar((C.transpose() * F).sum());


    return  EXIT_SUCCESS;
}
