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
    index_t npts = 100;

    gsCmdLine cmd("Computes patches from structured (tensor-product) data samples by fitting.");
    cmd.addPlainString("filename", "File containing multipatch input (.xml).", filename);
    cmd.addReal  ("t","tolerance","Tolerance for identifing patch interfaces", tol);
    cmd.addInt   ("d", "degree", "Degree of B-splines for reparameterization", degree);
    cmd.addInt   ("k", "knots", "Number of interior knots for reparameterization", nknots);
    cmd.addInt   ("N", "npts", "Number of points for sampling", npts);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<>::uPtr mp0 = gsReadFile<>(filename);
    mp0->computeTopology();
    gsInfo <<" Got"<< *mp0 <<" \n" ;


    // STEP 0: sample all patches to get linear patches.
    // (already the input is ASSUMED linear for now)
    // sample each patch on a grid
    // create the linear patches with the samples as coefficients + the topology of the initial one   
    gsMatrix<> ab, pts, eval;
    gsVector<> a, b;
    gsVector<unsigned> np;
    gsMultiPatch<> mp = *mp0;
    gsMultiPatch<> mp_par = *mp0;
    for (size_t p=0; p!=mp0->nPatches(); p++)
    {
        ab = mp0->patch(p).support();
        a = ab.col(0);
        b = ab.col(1);
        np = uniformSampleCount(a, b, npts);
        pts = gsPointGrid(a, b, np);

        mp0->patch(p).eval_into(pts,eval);

        // maybe also choosae knot-parametrization (uniform, centripetal or others...)
        gsKnotVector<> kv0(0, 1, np[0]-2, 2, 1, 1);
        gsKnotVector<> kv1(0, 1, np[1]-2, 2, 1, 1);
        gsTensorBSplineBasis<2,real_t> bbasis(kv0,kv1);

        mp.patch(p) = *bbasis.makeGeometry(eval.transpose()).release();
        mp_par.patch(p) = *bbasis.makeGeometry(pts.transpose()).release();
    }
    // mp.topology() = mp0->topology();

    gsDebugVar(mp.patch(0));
    gsWriteParaview(mp,"mp",1000,true);

    // todo: mp = mp0->approximateLinearly();
    
    // STEP 1: Get curve network with merged linear interfaces
    mp.constructInterfaceRep();
    mp.constructBoundaryRep();
    auto & irep = mp.interfaceRep();
    auto & brep = mp.boundaryRep();
    gsInfo <<" irep "<< irep.size() <<" \n" ;
    gsInfo <<" brep "<< brep.size() <<" \n" ;

    // outputing...
    gsMultiPatch<> crv_net;
    for (auto it = irep.begin(); it!=irep.end(); ++it)
        crv_net.addPatch((*it->second));
    for (auto it = brep.begin(); it!=brep.end(); ++it)
        crv_net.addPatch((*it->second));
    // gsFileData<> fd;
    // fd<< crv_net ;
    // fd.dump("linear_output");
    // gsInfo <<"Resulting file is written out to linear_output.xml\n" ;
    gsWriteParaview(crv_net,"crv_net",1000,true);
    //end outputing
    
    //STEP 3: Fit curve network with B-splines of degree \a d and \a k interior knots
    //parametrizePts
    gsKnotVector<> kv(0, 1, nknots, degree+1, 1, degree);
    gsDebugVar(kv);
    gsBSplineBasis<> fbasis(kv);
    gsFitting<> cfit;

    std::vector<gsMultiPatch<> > pbdr(mp.nPatches()); // patch boundary curves (4 for each)
    crv_net.clear();
    gsMatrix<> uv, xyz;
    // Can we make a nice initialization of the multipatch? Such that we can put the sides directly in the side index
    for (size_t p=0; p!=mp.nPatches(); p++)
        for (index_t k=0; k!=4; k++)
            pbdr[p].addPatch(mp.patch(0));


    index_t k=0;
    for (auto it = irep.begin(); it!=irep.end(); ++it, k++)
    {
        const gsBSpline<> & crv = static_cast<const gsBSpline<> &>(*it->second);
        // note: free to choose a parameterization, here we just take the knots
        // uv  = crv.knots().asMatrix().middleCols(1,crv.numCoefs());
        //
        // THIS IS VERY WRONG
        xyz = crv.coefs().transpose();
        xyz.conservativeResize(xyz.rows(),xyz.cols()-1);
        uv = gsVector<>::LinSpaced(xyz.cols(),0,1).transpose();

        gsDebugVar( xyz );
        gsDebugVar( uv );
        gsDebugVar(crv.knots().asMatrix());
        cfit = gsFitting<>(uv,xyz,fbasis);
        cfit.compute();
        pbdr[it->first.first() .patch].patch(it->first.first() .side()-1) = *cfit.result() ;
        pbdr[it->first.second().patch].patch(it->first.second().side()-1) = *cfit.result() ;
        crv_net.addPatch( *cfit.result() );

        if (k==1)
        {
            gsWriteParaview(*cfit.result(),"result");
            gsMatrix<> coefs = cfit.result()->coefs().transpose();
            gsWriteParaviewPoints(coefs,"crv_fit_pts");
            gsWriteParaviewPoints(xyz,"original_pts");
        }
    }
    for (auto it = brep.begin(); it!=brep.end(); ++it)
    {
        const gsBSpline<> & crv = static_cast<const gsBSpline<> &>(*it->second);
        // note: free to choose a parameterization, here we just take the knots
        uv  = crv.knots().asMatrix().middleCols(1,crv.numCoefs());
        xyz = crv.coefs().transpose();
        cfit = gsFitting<>(uv,xyz,fbasis);
        cfit.compute();
        pbdr[it->first.patch].patch(it->first.side()-1) = *cfit.result() ;
        crv_net.addPatch( *cfit.result() );
    }

    gsWriteParaview(crv_net.patch(1),"crv_fit",1000,false,true);
    gsMatrix<> coefs = crv_net.patch(1).coefs().transpose();
    gsWriteParaviewPoints(coefs,"crv_fit_pts");

    // fd.clear();
    // fd<< crv_net ;
    // fd.dump("curve_output");
    // gsInfo <<"Resulting file is written out to curve_output.xml\n" ;

    gsTensorBSplineBasis<2> sbasis(kv,kv);
    gsFitting<> sfit;
    gsMultiPatch<> mp_res;
    //STEP 4: fit interior points of each patch with boundary constraints being the curves..
    for (size_t p=0; p!=pbdr.size(); p++)
    {
        uv  = mp_par.patch(p).coefs().transpose();
        xyz = mp.patch(p).coefs().transpose();
        sfit = gsFitting<>(uv,xyz,sbasis);

        std::vector<gsGeometry<real_t> * > prescribedCurves;
        std::vector<boxSide> prescribedSides;

        for (index_t s=0; s!=4; s++) //sides
        {
            prescribedSides.push_back(boxSide(s+1));
            prescribedCurves.push_back(&pbdr[p].patch(s));
        }
        sfit.setConstraints(prescribedSides, prescribedCurves);

        sfit.compute();
        gsWriteParaview(*sfit.result(),"fit");
        mp_res.addPatch(*sfit.result());
    }
   
    gsWriteParaview(mp_res,"final",1000,true);

    return EXIT_SUCCESS;
}
