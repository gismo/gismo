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
    fn = "../filedata/surfaces/simple.xml";
    index_t addRef = 1;

    gsCmdLine cmd("Input: .xml file with gsTensorProductBspline geometry.");
    cmd.addInt("r", "ref", "Uniform refinement", addRef);
    cmd.addString("m", "measuredata", "Input measure data", fn);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (fn.empty() )
    {
        gsInfo<< cmd.getMessage();
        gsInfo<<"\nType "<< argv[0]<< " -h, to get the list of command line options.\n";
        return 0;
    }

    gsFileData<>  data( fn );
    data.save();

    gsInfo<< "Measurement file contains "<< data.numTags() <<" Geometric objects.\n";
    memory::unique_ptr<gsGeometry<> > surf = data.getFirst< gsGeometry<> >();
    gsInfo << *surf << "\n";
    
    gsWriteParaview(*surf, "surf", 10000, false, true);

    gsMultiBasis<> mb(surf->basis()); // get the basis

    // Method 1.
    gsExprEvaluator<> ev1;
    ev1.setIntegrationElements(mb);

    typedef typename gsExprEvaluator<>::geometryMap geometryMap;
    geometryMap Gmeasure1 = ev1.getMap(*surf);

    gsVector<> ez_vec(3);
    ez_vec<<0,0,1;
    gsConstantFunction<> ez_fun(ez_vec,2);
    auto ez = ev1.getVariable(ez_fun);

    real_t volume1 = ev1.integral(Gmeasure1.tr() * ez * meas(Gmeasure1));

    gsInfo<< "Volume_1 = " << volume1 << "\n\n";
    gsInfo<< "---------------------\n\n";


    // Method 2.
    //1. Take the projection in (x,y). This will be a planar B-spline, lets call it G.
    //2. Take a spline function only being the z-part. This is a scalar B-spline, let's call it f.
    //3. Use this planar domain as your geometry map in the gsExprEvaluator
    //4. Compute the integral:    integral( f * meas(G) ) .
    //   To get better accuracy for the computed integral, you can uniformly refine the G if needed

    real_t volume2 = 0.;

    // 1. Take the projection in (x,y). This will be a planar B-spline, lets call it G.
    gsMatrix<> dcoefs(surf->coefs().rows(), 2);
    dcoefs.col(0) = surf->coefs().col(0);
    dcoefs.col(1) = surf->coefs().col(1);

    // assemble integration domain
    auto domainG = surf->basis().makeGeometry(dcoefs);

    //2. Take a spline function only being the z-part. This is a scalar B-spline, let's call it f.
    gsMatrix<> zcoefs(surf->coefs().rows(), 1);
    zcoefs = surf->coefs().col(2);

    // define integrand
    auto f = surf->basis().makeGeometry(zcoefs);


    gsMatrix<> domainCoefs = domainG->coefs().transpose();
    gsWriteParaviewPoints(domainCoefs, "domainCoefs");
    gsWriteParaview(*domainG, "domain");
            
    
    //3. Use this planar domain as your geometry map in the gsExprEvaluator
    gsExprEvaluator<> ev;
    ev.setIntegrationElements(mb);
    typedef typename gsExprEvaluator<>::geometryMap geometryMap;
    geometryMap Gmeasure = ev.getMap(*domainG);

    //4. Compute the integral:    integral( f * meas(G) ) .
    auto fz = ev.getVariable(*f);
    volume2 = ev.integral(fz * meas(Gmeasure));
    gsInfo<< "Volume_2 = " << volume2 << "\n";


     // If we consider the ingral to be always the [0,1]^2 domain;
    real_t volume3 = ev.integral(fz);
    gsInfo<< "Volume_3 = " << volume2 << "\n";

    return EXIT_SUCCESS;
}
