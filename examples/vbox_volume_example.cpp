/** @file vbox_volume_example.cpp
    @brief Testing file reading and writing
    This file is part of the G+Smo library.
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    Author(s): S. Imperatore


    1. Take the projection in (x,y). This will be a planar B-spline, lets call it G.
    2. Take a spline function only being the z-part. This is a scalar B-spline, let's call it f.
    3. Use this planar domain as your geometry map in the gsExprEvaluator
    4. Compute the integral:    integral( f * meas(G) ) .
        To get better accuracy for the computed integral, you can uniformly refine the G if needed
*/

#include <iostream>

#include <gismo.h>
using namespace gismo;

int main(int argc, char *argv[])
{
    std::string fn;
    fn = "../filedata/vboxdata/posB_lpsp_fit.xml";

    gsCmdLine cmd("Input: .xml file with gsTensorProductBspline geometry.");
    cmd.addString("m", "measuredata", "Input measure data", fn);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (fn.empty() )
    {
        gsInfo<< cmd.getMessage();
        gsInfo<<"\nType "<< argv[0]<< " -h, to get the list of command line options.\n";
        return 0;
    }


    gsMatrix<> pltCoefs;

    gsFileData<>  data( fn );
    data.save();


    gsInfo<< "Measurement file contains "<< data.numTags() <<" gsTensorBSpline objects.\n";
    memory::unique_ptr<gsGeometry<> > surf = data.getFirst< gsGeometry<> >();
    
    // resetting the coordinate system for the integral computation;
    surf->coefs().col(2).array() *= -1;
    surf->coefs().col(2).array() += 1.; // to be scaled properly, accordin to the the data magnitude

    surf->coefs().col(0).array() *= -1;
    surf->coefs().col(0).array() += 1; // to be scaled properly, accordin to the the data magnitude
    
    gsInfo << *surf << "\n";
    gsWriteParaview(*surf, "surf", 10000, false, true);

    gsFileData<> tmp;
    tmp << *surf;
    tmp.dump("posB");

    pltCoefs = surf->coefs().transpose();
    gsWriteParaviewPoints(pltCoefs, "surf_net");

    gsMultiBasis<> mb(*surf);
    gsExprEvaluator<> ev;
    ev.setIntegrationElements(mb);

    typedef typename gsExprEvaluator<>::geometryMap geometryMap;
    geometryMap Gmeasure = ev.getMap(*surf);

    gsVector<> ez_vec(3);
    ez_vec<<0,0,1;
    gsConstantFunction<> ez_fun(ez_vec,2);
    auto ez = ev.getVariable(ez_fun);

    real_t volume = ev.integral(Gmeasure.tr() * ez * meas(Gmeasure));
    // real_t volume = ev.integral(Gmeasure.norm()*meas(Gmeasure));

    gsInfo<< "Volume = " << volume << "\n";

    return EXIT_SUCCESS;
}
