/** @file L2projection_example.cpp

    @brief This file contains an example of using the L2Projection class in the G+Smo library to project a geometry and a function onto a finer basis.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
**/

/* 
    TO DO (HMV):
    * Make unit test
    * Add multipatch case
    * Add case with different integration basis
 */

#include <gismo.h>

using namespace gismo;

int main(int argc, char* argv[])
{
    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    // cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    // cmd.addString( "f", "file", "Input XML file", fn );
    // cmd.addSwitch("binary", "Use B64 encoding for Paraview", export_b64);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    // Make a geometry
    gsMultiPatch<> mp;
    mp.addPatch(gsNurbsCreator<>::NurbsAnnulus());

    // Make a basis
    gsMultiBasis<> mb(mp);
    mb.uniformRefine();

    gsMultiBasis<> mb_coarse(mp);

    // Coefficients object
    gsMatrix<> coefs_SP, coefs_MP;

    // Error
    real_t error;

    gsInfo<<"---------------------------------------------------------------------------------------\n";
    gsInfo<<"-------------------------------Geometry projection-------------------------------------\n";
    gsInfo<<"---------------------------------------------------------------------------------------\n";

    // Project the geometry on the finer basis
    // Single patch version
    gsInfo<<"Projection 1, error = "<<gsL2Projection<real_t>::project(mb.basis(0),mp.patch(0),coefs_SP)<<"\n";
    // Multipatch version
    gsInfo<<"Projection 2, error = "<<gsL2Projection<real_t>::project(mb,mp,coefs_MP)<<"\n";
    // Check the difference
    GISMO_ENSURE((coefs_SP-coefs_MP).norm() < 1e-10,"Difference between single patch and multipatch projection is too large");

    // Construct the geometry
    gsGeometry<>::uPtr geom = mb.basis(0).makeGeometry(give(coefs_SP.reshape(mb.basis(0).size(),mp.targetDim())));

    // Compute the error
    gsExprEvaluator<> ev;
    ev.setIntegrationElements(mb);
    auto G = ev.getMap(mp);
    auto ori = ev.getVariable(mp);
    auto proj= ev.getVariable(*geom);
    error = ev.integral((ori-proj).sqNorm()*meas(G));
    gsDebugVar(error);


    gsInfo<<"---------------------------------------------------------------------------------------\n";
    gsInfo<<"-------------------------------Function projection-------------------------------------\n";
    gsInfo<<"---------------------------------------------------------------------------------------\n";

    std::vector<std::string> expressions = {"x^2*y","y^2*x","sin(x)*sin(y)","x","y"};
    gsFunctionExpr<> f(expressions,2);

    // Project the geometry on the finer basis
    // Single patch version
    gsInfo<<"Projection 1, error = "<<gsL2Projection<real_t>::project(mb.basis(0),mp.patch(0),f,coefs_SP)<<"\n";
    // Multipatch version
    gsInfo<<"Projection 2, error = "<<gsL2Projection<real_t>::project(mb,mp,f,coefs_MP)<<"\n";
    // Check the difference
    GISMO_ENSURE((coefs_SP-coefs_MP).norm() < 1e-10,"Difference between single patch and multipatch projection is too large");

    gsGeometry<>::uPtr ff = mb.basis(0).makeGeometry(give(coefs_SP.reshape(mb.basis(0).size(),f.targetDim())));

    auto fori = ev.getVariable(f);
    auto fproj= ev.getVariable(*ff);
    error = ev.integral((fori-fproj).sqNorm()*meas(G));
    gsDebugVar(error);


    // VECTOR


    return EXIT_SUCCESS;
}
