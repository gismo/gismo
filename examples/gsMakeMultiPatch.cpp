/** @file gsMakeMultiPatch.cpp

    @brief Computes correctly the boundaries and interfaces of
    a multipatch structure

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
    std::string filename("planar/lshape2d_3patches_tens.xml");
    real_t tol = 1e-5;
    real_t gtol = 1e-6;
    bool reparam = false, gaps = true;
    index_t method = 0, nknots = 5, degree = 3;

    gsCmdLine cmd("Computes the topology of a set of patches, identifing interfaces and boundaries.");
    cmd.addPlainString("filename", "File containing multipatch input (.xml).", filename);
    cmd.addReal  ("t","tolerance","Tolerance for identifing patch interfaces", tol);
    cmd.addReal  ("g","gap-tolerance","Tolerance for closing gaps", gtol);
    cmd.addInt   ("d", "degree", "Degree of B-splines for reparameterization", degree);
    cmd.addInt   ("k", "knots", "Number of interior knots for reparameterization", nknots);
    cmd.addSwitch("reparam", "Reparameterize all patches using a fixed degree and number of knots", reparam);
    cmd.addSwitch("nogaps", "Close any gaps along interfaces upto tolerance \'gtol\'", gaps);

    cmd.addInt   ("m", "method", "method", method);
        
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<>::uPtr mp = gsReadFile<>(filename);
    gsInfo <<" Got"<< *mp <<" \n" ;

    size_t d = mp->parDim();
    for ( size_t i = 0; i< mp->nPatches(); ++i )
    {
        gsTensorBSplineBasis<2> * bb =
            dynamic_cast<gsTensorBSplineBasis<2> *>(&mp->patch(i).basis().source());
        GISMO_ENSURE(nullptr!=bb, "Conversion error.");
        for ( size_t k = 0; k!=d; ++k )
            bb->component(k).knots().transform(0,1);
    }

    if (reparam) // Reparameterize?
    {
        gsKnotVector<> kv(0,1,nknots,degree+1,1);
        gsTensorBSplineBasis<2> bs(kv,kv);
        gsMultiPatch<> newmp;
        // Using interpolation
        if (0==method)
        {
            gsMatrix<> values, anchors = bs.anchors();
            for ( size_t i = 0; i< mp->nPatches(); ++i )
            {
                values = mp->patch(i).eval(anchors);
                newmp.addPatch( bs.interpolateData(values,anchors) );
            }
        }

        // Using Fitting

        if (1==method)
        {
            gsVector<index_t,2> sz; //str
            gsVector<unsigned> sz2(2);
            
            gsVector<> a(2); a.setZero();
            gsVector<> b(2); b.setOnes();
            gsMatrix<> param;
            
            for ( size_t i = 0; i< mp->nPatches(); ++i )
            {
                gsTensorBSplineBasis<2> * bb =
                    dynamic_cast<gsTensorBSplineBasis<2> *>(&mp->patch(i).basis().source());
                bb->size_cwise(sz);
                sz2[0]= sz[0];
                sz2[1]= sz[1];
                param = gsPointGrid(a,b,sz2);
                gsFitting<>  fitting(param, mp->patch(i).coefs().transpose(), bs);
                fitting.compute(0.0000001);
                newmp.addPatch( *fitting.result() );

                fitting.computeMaxNormErrors();
                gsInfo<< "Min error: "<< fitting.minPointError() <<"\n";
                gsInfo<< "Max error: "<< fitting.maxPointError() <<"\n";
            }
        }
        
        if (1==method)
        {
            gsVector<index_t,2> sz; //str
            gsVector<unsigned> sz2(2);
            
            gsVector<> a(2); a.setZero();
            gsVector<> b(2); b.setOnes();
            gsMatrix<> param;
            
            for ( size_t i = 0; i< mp->nPatches(); ++i )
            {
                gsTensorBSplineBasis<2> * bb =
                    dynamic_cast<gsTensorBSplineBasis<2> *>(&mp->patch(i).basis().source());
                bb->size_cwise(sz);
                sz2[0]= sz[0];
                sz2[1]= sz[1];
                param = gsPointGrid(a,b,sz2);
                gsFitting<>  fitting(param, mp->patch(i).coefs().transpose(), bs);
                fitting.compute(0.0000001);
                newmp.addPatch( *fitting.result() );

                fitting.computeMaxNormErrors();
                gsInfo<< "Min error: "<< fitting.minPointError() <<"\n";
                gsInfo<< "Max error: "<< fitting.maxPointError() <<"\n";
            }
        }

        mp->swap(newmp);
    }

    gsInfo <<"Computing"<< (gaps?" corner-based ":" ")<<"topology with tolerance = "<<tol<<"... \n" ;
    mp->computeTopology(tol, gaps);
    gsInfo << * mp <<"\n" ;
    //gsInfo << mp->detail() <<"\n";
    if (gaps) //Close gaps?
    {
        gsInfo <<"Closing gaps with tolerance = "<<gtol<<"... \n" ;
        mp->closeGaps(gtol);
        gsInfo <<"Computing topology with tolerance = "<<tol<<"... \n" ;
        mp->computeTopology(tol, false);
        mp->closeGaps(gtol);//close again for any previously missed interfaces
        gsInfo << * mp <<"\n" ;
        //gsInfo << mp->detail() <<"\n";
    }

    gsFileData<> fd;
    fd<< *mp ;
    fd.dump("makeMultipatch_output");    
    gsInfo <<"Resulting file is written out to makeMultipatch_output.xml\n" ;

    return EXIT_SUCCESS;
}
