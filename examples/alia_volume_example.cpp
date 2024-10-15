/** @file alia_volume_example.cpp
    @brief Testing file reading and writing
    This file is part of the G+Smo library.
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    Author(s): S. Imperatore
*/

#include <iostream>

#include <gismo.h>
using namespace gismo;

int main(int argc, char *argv[])
{
    std::string fc, fn;
    fc = "";
    fn = "../filedata/alia/posB_lpsp_fit.xml";

    gsCmdLine cmd("Input: .xml file with gsTensorProductBspline approximating the trash-bin content.");
    cmd.addString("m", "measuredata", "Input measure data", fn);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (fn.empty() )
    {
        gsInfo<< cmd.getMessage();
        gsInfo<<"\nType "<< argv[0]<< " -h, to get the list of command line options.\n";
        return 0;
    }


    gsMatrix<> pltCoefs;

    if (fc.length() > 0)
    {
      gsFileData<>  data_calib( fc );
      data_calib.save();

      gsInfo<< "Calibration file contains "<< data_calib.numTags() <<" gsTensorBSpline objects.\n";
      memory::unique_ptr<gsTensorBSpline<2, real_t> > calib = data_calib.getFirst< gsTensorBSpline<2, real_t> >();
      gsInfo << *calib << "\n";
      gsWriteParaview(*calib, "calib", 10000, false, true);
      pltCoefs = calib->coefs().transpose();
      gsWriteParaviewPoints(pltCoefs, "calib_net");
    }

    gsFileData<>  data( fn );
    data.save();


    gsInfo<< "Measurement file contains "<< data.numTags() <<" gsTensorBSpline objects.\n";
    memory::unique_ptr<gsTHBSpline<2, real_t> > measure = data.getFirst< gsTHBSpline<2, real_t> >();
    
    // resetting the coordinate system for the integral computation;
    measure->coefs().col(2).array() *= -1;
    measure->coefs().col(2).array() += 1.;

    measure->coefs().col(0).array() *= -1;
    measure->coefs().col(0).array() += 1;
    
    gsInfo << *measure << "\n";
    gsWriteParaview(*measure, "measure", 10000, false, true);

    gsFileData<> tmp;
    tmp << *measure;
    tmp.dump("posB");

    pltCoefs = measure->coefs().transpose();
    gsWriteParaviewPoints(pltCoefs, "measure_net");


    // Assuming the bases of both geometries are the same
    gsMultiBasis<> mb(*measure);
    gsExprEvaluator<> ev;
    ev.setIntegrationElements(mb);

    typedef typename gsExprEvaluator<>::geometryMap geometryMap;
    geometryMap Gmeasure = ev.getMap(*measure);

    gsVector<> ez_vec(3);
    ez_vec<<0,0,1;
    gsConstantFunction<> ez_fun(ez_vec,2);
    auto ez = ev.getVariable(ez_fun);

    real_t volume = ev.integral(Gmeasure.tr() * ez * meas(Gmeasure));
    // real_t volume = ev.integral(Gmeasure.norm()*meas(Gmeasure));


    gsInfo<<volume<<"\n";

    return EXIT_SUCCESS;
}
