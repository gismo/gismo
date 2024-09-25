/** @file lofting_exm_pipeline.cpp

    @brief Demonstrates adaptive fitting of data samples

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore
*/

#include <gismo.h>
#include <gsModeling/gsLofting.h>

using namespace gismo;

gsMultiPatch<> getMPfromGeom(gsFileData<> & data);

template<class T>
gsMatrix<T> parameterize_points1D(const gsMatrix<T> & xyz, T alpha = 0.5)
{
    index_t N = xyz.cols();
    gsMatrix<T> uv(1,N);

    uv(0,0)   = 0;
    uv(0,N-1) = 1;
    std::vector<T> distances(N-1);
    for (index_t k=1; k!=N; k++)
        distances.at(k-1) = math::pow((xyz.col(k)-xyz.col(k-1)).norm(),alpha);
    T d = std::accumulate(distances.begin(),distances.end(), T(0));
    for (index_t k=1; k!=N-1; k++)
        uv(0,k) = uv(0,k-1) + distances.at(k) / d;

    return uv;
}




int main(int argc, char *argv[])
{
    // Options with default values
    
    std::string fn = "../filedata/fitting/vee_rail6_pts.xml";

    real_t lambda = 1e-07;
    index_t numSamples = 1000;
    index_t deg_v = 2;
    index_t deg_u = 2;
    index_t num_u_knots = 5;
    real_t cm = 1.;

    // Reading options from the command line
    gsCmdLine cmd("Loft parametrized B-spline cuves. Expected input file is an XML "
            "file containing the data points to be fitted by curves and then lofted.");
    
    cmd.addString("f", "filedata", "Input sample data", fn);
    cmd.addReal("l", "lambda", "smoothing coefficient", lambda);
    cmd.addInt("s", "samples", "number of plotting samples", numSamples);
    cmd.addInt("u", "udeg", "degree in u direction", deg_u);
    cmd.addInt("v", "vdeg", "degree in v direction", deg_v);
    cmd.addInt("", "nu", "number of knots in u direction", num_u_knots);
    cmd.addReal("m", "method", "curves parameterization method: (0) unfirom; (0.5) centripetal; (1) chord-length", cm);
    

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<real_t> container;

    gsKnotVector<> ku (0., 1., num_u_knots, deg_u+1) ;

    gsFileData<> fd(fn);
    gsMatrix<> xyz, uv;
    

    gsInfo << "Numdata: " << fd.numData() << "\n";


    std::vector< memory::unique_ptr< gsMatrix<> > >  mContainer = fd.getAll< gsMatrix<> >(); 

    for(index_t j=0; j < mContainer.size(); j ++)
    {
        fd.getId<gsMatrix<> >(j+1, xyz );

        gsInfo << j << "j-th input data:\n";
        gsInfo << xyz.rows() << " x " << xyz.cols() << "\n";

        uv = parameterize_points1D(xyz, cm);

        gsInfo << uv.rows() << " x " << uv.cols() << "\n";

        //xyz.transposeInPlace();
        //uv.transposeInPlace();
        gsCurveFitting<real_t> fitter( uv, xyz, ku);
        fitter.compute(lambda);

        gsInfo << fitter.curve() << "\n";
        container.addPatch(fitter.curve());
    }
    

    gsInfo << container << "\n";
    gsInfo << "number of lofting curves: " << container.nPatches() << "\n";

    gsWriteParaview(container, "section_curves", numSamples);
    
    // choose the degree, parameters and knot-vector in the v-direction.
    // degree must be less than #section curves

    index_t num_v_knots = container.nPatches() - (deg_v + 1);
    gsKnotVector<> kv (0., 1., num_v_knots, deg_v+1) ;


    // Create lofting object
    // gsLofting<real_t> loft(container);
    // gsLofting<real_t> loft(container, deg_v);
    gsLofting<real_t> loft(container, kv);
    loft.compute();

    // This is for outputing an XML file, if requested
    gsFileData<> out;
    out << *loft.result() ;

    out.dump("lofting_out");


    return 0;
}


gsMultiPatch<> getMPfromGeom(gsFileData<> & data)
{
    gsMultiPatch<> container;
    std::vector< memory::unique_ptr< gsGeometry<> > >  geoContainer = data.getAll< gsGeometry<> >();

    for(index_t el = 0; el < geoContainer.size(); el++)
        container.addPatch(*geoContainer[el]);

    return container;
}
