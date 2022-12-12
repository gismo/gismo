/*  @file fitting_msplines.cpp

    @brief Implementation of point cloud data fitting with
           multipatch splines functions

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: M. Marsala

*/

#include <iostream>
#include <time.h>
#include <gismo.h>


using namespace gismo;


int main(int argc, char* argv[])
{
    bool plot = false;
    index_t maxPcIter = 1;
    std::string filename = "fitting/chess.xml";

    int np = 700;
    real_t lambda = 0;
    std::string fn = "fitting/acc_chess_ptcloud.xml";

    gsCmdLine cmd("Fit parametrized sample data with a surface patch. Expected input file is an XML "
        "file containing two matrices (<Matrix>), with \nMatrix id 0 : contains a 2 x N matrix. "
        "Every column represents a (u,v) parametric coordinate\nMatrix id 1 : contains a "
        "3 x N matrix. Every column represents a point (x,y,z) in space.");
    cmd.addString("g", "filename", "File containing mp geometry and basis (.xml).", filename);
    cmd.addString("d", "data", "Input sample data", fn);
    cmd.addReal("l", "lambda", "smoothing parameter", lambda);
    cmd.addInt("c", "parcor", "Steps of parameter correction", maxPcIter);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsFileData<> fd_in(fn);
    gsMatrix<> Mpar, fval;
    gsMatrix<index_t> offset, pId;
    fd_in.getId<gsMatrix<> >(0, Mpar);
    fd_in.getId<gsMatrix<> >(1, fval);
    if (fd_in.hasId(2))
    {
        fd_in.getId<gsMatrix<index_t> >(2, offset);
    }
    if (fd_in.hasId(3))
    {
        fd_in.getId<gsMatrix<index_t> >(3, pId);

    }

    

    // gsVector<index_t>(offset);
    
    // The file data
    gsFileData<> data( filename );
    gsFunctionExpr<> f;
    gsMultiPatch<> mp;
    gsMultiBasis<> mb;
    gsSparseMatrix<> cf;
    gsMappedBasis<2, real_t> mbasis;

    // Load data as multipatch structure

    data.getFirst(mp);
    data.getFirst(mb);
    data.getFirst(cf);

    mbasis.init(mb, cf);

    gsInfo << "Number of patches: " << mp.nPatches() << "\n";
    //gsInfo << "Number of points per patch: " << nd << "\n";
    gsInfo << "Number of basis functions: " << mbasis.size() << "\n";
    //gsInfo << "Basis" << mb << "\n";

    //create the fitting object
    gsInfo << "//////////////////////////////////////////////////////////" << "\n";
    gsInfo << "///////   Creating the multipatch fitting object   ///////" << "\n";
    gsInfo << "//////////////////////////////////////////////////////////" << "\n";


    gsFitting<> fitting(Mpar, fval, offset, mbasis);


    gsInfo << "Fit class created" << "\n";

    fitting.compute(lambda); // smoothing parameter lambda

    fitting.parameterCorrection(maxPcIter);

    gsInfo << "I computed the fitting" << "\n";

    gsMappedSpline<2, real_t>  test;
    test = fitting.mresult();

    std::vector<real_t> errors;
    fitting.get_Error(errors, 0);

    real_t min_error = *std::min_element(errors.begin(), errors.end());
    real_t max_error = *std::max_element(errors.begin(), errors.end());

    gsInfo << "Min error: " << min_error << "\n";
    gsInfo << "L_inf: " << max_error << "\n";

    real_t error;

    fitting.computeApproxError(error, 0);

    gsInfo << "L_2: " << error << "\n";







    gsFileData<> newdata;
    gsMultiPatch<> surf = test.exportToPatches();
    newdata << surf;
    newdata.save("fitting_result");

    if (plot)
    {
        gsInfo << "Plotting in Paraview..." << "\n";
        gsWriteParaviewPoints(fval, "point_data");
        gsWriteParaview(surf, "multipatch_spline", np);
        gsFileManager::open("multipatch_spline.pvd");
    }

    return EXIT_SUCCESS;
}
