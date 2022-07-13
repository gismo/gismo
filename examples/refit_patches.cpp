/** @file gsFitPatches.cpp

    @brief Computes patches from structured (tensor-product) data samples by fitting.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    std::string filename("domain2d/yeti_mp2.xml");
    real_t tol = 1e-5;
    index_t nknots = 5, degree = 3;

    gsCmdLine cmd("Computes patches from structured (tensor-product) data samples by fitting.");
    cmd.addPlainString("filename", "File containing multipatch input (.xml).", filename);
    cmd.addReal  ("t","tolerance","Tolerance for identifing patch interfaces", tol);
    cmd.addInt   ("d", "degree", "Degree of B-splines for reparameterization", degree);
    cmd.addInt   ("k", "knots", "Number of interior knots for reparameterization", nknots);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<>::uPtr mp0 = gsReadFile<>(filename);
    gsInfo <<" Got"<< *mp0 <<" \n" ;


    // STEP 0: sample all patches to get linear patches.
    // (already the input is ASSUMED linear for now)
    // sample each patch on a grid
    // create the linear patches with the samples as coefficients + the topology of the initial one   
    gsMultiPatch<>::Ptr mp = mp0->clone();//sampling...
    // todo: mp = mp0->approximateLinearly();
    //maybe also choosae knot-parametrization (uniform, centripetal or others...)
    
    // STEP 1: Get curve network with merged linear interfaces
    mp->computeTopology();
    mp->constructInterfaceRep();
    mp->constructBoundaryRep();
    auto & irep = mp->interfaceRep();
    auto & brep = mp->boundaryRep();
    gsInfo <<" irep "<< irep.size() <<" \n" ;
    gsInfo <<" brep "<< brep.size() <<" \n" ;

    // outputing...
    gsMultiPatch<> crv_net;
    for (auto it = irep.begin(); it!=irep.end(); ++it)
        crv_net.addPatch((*it->second));
    for (auto it = brep.begin(); it!=brep.end(); ++it)
        crv_net.addPatch((*it->second));
    gsFileData<> fd;
    fd<< crv_net ;
    fd.dump("linear_output");
    gsInfo <<"Resulting file is written out to linear_output.xml\n" ;
    //end outputing
    
    //STEP 3: Fit curve network with B-splines of degree \a d and \a k interior knots
    //parametrizePts
    gsKnotVector<> kv(0, 1, nknots, degree+1, 1, degree);
    gsBSplineBasis<> fbasis(kv);
    gsFitting<> cfit;

    std::vector<gsMultiPatch<> > pbdr(mp->nPatches()); // patch boundary curves (4 for each)
    crv_net.clear();
    gsMatrix<> uv, xyz;
    for (auto it = irep.begin(); it!=irep.end(); ++it)
    {
        const gsBSpline<> & crv = static_cast<const gsBSpline<> &>(*it->second);
        // note: free to choose a parameterization, here we just take the knots
        uv  = crv.knots().asMatrix().middleCols(1,crv.numCoefs());
        xyz = crv.coefs().transpose();
        gsDebugVar( xyz );
        gsDebugVar( uv );
        cfit = gsFitting<>(uv,xyz,fbasis);
        cfit.compute();
        pbdr[it->first.first() .patch].addPatch( *cfit.result() );
        pbdr[it->first.second().patch].addPatch( *cfit.result() );
        crv_net.addPatch( *cfit.result() );
    }
    for (auto it = brep.begin(); it!=brep.end(); ++it)
    {
        const gsBSpline<> & crv = static_cast<const gsBSpline<> &>(*it->second);
        // note: free to choose a parameterization, here we just take the knots
        uv  = crv.knots().asMatrix().middleCols(1,crv.numCoefs());
        xyz = crv.coefs().transpose();
        cfit = gsFitting<>(uv,xyz,fbasis);
        cfit.compute();
        pbdr[it->first.patch].addPatch( *cfit.result() );
        crv_net.addPatch( *cfit.result() );
    }

    fd.clear();
    fd<< crv_net ;
    fd.dump("curve_output");
    gsInfo <<"Resulting file is written out to curve_output.xml\n" ;


    //STEP 4: fit interior points of each patch with boundary constraints being the curves..

   
    return EXIT_SUCCESS;
}
