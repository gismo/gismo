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

#ifdef gsOpennurbs_ENABLED
#include "gsOpennurbs/gsWriteOpenNurbs.h"
#endif


using namespace gismo;



gsSparseMatrix<> makeC0basis(const gsMultiBasis<>& mb)
{
    gsDofMapper m_dofMapper;
    mb.getMapper(true, m_dofMapper);
    index_t m_globals = m_dofMapper.size();
    index_t m_locals = m_dofMapper.mapSize();
    gsSparseMatrix<> matrix(m_locals, m_globals);
    index_t local = 0;
    for (size_t patch = 0; patch < m_dofMapper.numPatches(); ++patch)
    {
        const size_t patchSize = (patch != m_dofMapper.numPatches() - 1) ?
            m_dofMapper.offset(patch + 1) - m_dofMapper.offset(patch) :
            m_locals - m_dofMapper.offset(m_dofMapper.numPatches() - 1);
        for (size_t i = 0; i < patchSize; ++i)
        {
            const index_t global = m_dofMapper.index(i, patch);
            matrix.at(local, global) = 1;
            local++;
        }
    }
    matrix.makeCompressed();
    return matrix;
}


int main(int argc, char* argv[])
{
    bool bernstein = false,  plot = false;
    index_t maxPcIter = 1;
    std::string filename = "fitting/chess.xml";

    int np = 7000;
    real_t scale_par = 1;
    real_t lambda = 0;
    std::string fn = "fitting/acc_chess_ptcloud.xml";

    gsCmdLine cmd("Fit parametrized sample data with a surface patch. Expected input file is an XML "
        "file containing two matrices (<Matrix>), with \nMatrix id 0 : contains a 2 x N matrix. "
        "Every column represents a (u,v) parametric coordinate\nMatrix id 1 : contains a "
        "3 x N matrix. Every column represents a point (x,y,z) in space.");
    cmd.addString("g", "filename", "File containing mp geometry and basis (.xml).", filename);
    cmd.addString("d", "data", "Input sample data", fn);
    cmd.addReal("s", "scale_par", "Scale parameter for ptscloud", scale_par);
    cmd.addReal("l", "lambda", "smoothing parameter", lambda);
    cmd.addInt("c", "parcor", "Steps of parameter correction", maxPcIter);
    cmd.addSwitch("bb", "Use Bernstein basis", bernstein);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsFileData<> fd_in(fn);
    gsMatrix<> Mpar, fval;
    gsMatrix<index_t> offset, pId;
    fd_in.getId<gsMatrix<> >(0, Mpar);
    fd_in.getId<gsMatrix<> >(1, fval);

   


    if (fval.cols() == 3)
    {
        Mpar = Mpar.transpose();
        fval = fval.transpose();
    }

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
    typedef gsExprAssembler<>::space       space;

    // Load data as multipatch structure

    data.getFirst(mp);
    data.getFirst(mb);
    

    //gsMappedBasis<2, real_t> bb2;
    if (bernstein) { 
        cf = makeC0basis(mb);
    }
    else {
        data.getFirst(cf);
    }

    mbasis.init(mb, cf);


    gsInfo << "Number of patches: " << mp.nPatches() << "\n";
    //gsInfo << "Number of points per patch: " << nd << "\n";
    gsInfo << "Number of basis functions: " << mbasis.size() << "\n";
    //gsInfo << "Basis" << mb << "\n";

    //create the fitting object
    gsInfo << "//////////////////////////////////////////////////////////" << "\n";
    gsInfo << "///////   Creating the multipatch fitting object   ///////" << "\n";
    gsInfo << "//////////////////////////////////////////////////////////" << "\n";

    fval = scale_par * fval;

    gsStopwatch timer;

    timer.restart();

    gsFitting<> fitting(Mpar, fval, offset, mbasis);


    gsInfo << "Fit class created" << "\n";

    fitting.compute(lambda); // smoothing parameter lambda

    fitting.parameterCorrection(1e-12,maxPcIter);

    timer.stop();

    gsInfo << "Elapsed time for the solving: " << timer.elapsed() << " sec" << "\n";

    gsInfo << "I computed the fitting" << "\n";
   


    gsMappedSpline<2, real_t>  test;
    test = fitting.mresult();

    std::vector<real_t> errors;
    //fitting.get_Error(errors, 0);

    fitting.computeErrors();
    errors = fitting.pointWiseErrors();

    real_t ave_linf;
    ave_linf = 0;

   /* for (index_t j = 0; j < fval.cols(); j++) {

        ave_linf += errors[j];

    }*/

    gsMatrix<>ptswithcolors(4, fval.cols());

    ptswithcolors.topRows(3) = fval;
    ptswithcolors.row(3) = gsAsConstMatrix<>(errors);
    

    real_t min_error = *std::min_element(errors.begin(), errors.end());
    real_t max_error = *std::max_element(errors.begin(), errors.end());



    gsInfo << "Min error: " << min_error << "\n";
    gsInfo << "L_inf: " << max_error << "\n";
    //gsInfo << "Ave L_inf: " << ave_linf / fval.cols() << "\n";

    real_t error;

    fitting.computeApproxError(error, 0);

    gsInfo << "RMSE: " << std::sqrt(error/fval.cols()) << "\n";  








    gsFileData<> newdata;
    gsMultiPatch<> surf = test.exportToPatches();
    gsInfo << "Write back to mp.3dm\n";
    extensions::writeON_MultiPatch(surf,"mp");
    newdata << surf;
    newdata << mb;
    newdata << cf;
    newdata.save("fitting_result");
    gsWriteParaviewPoints(ptswithcolors, "point_data");

    if (plot)
    {
        gsInfo << "Plotting in Paraview..." << "\n";
        gsWriteParaviewPoints(ptswithcolors, "point_data");
        gsWriteParaview(surf, "multipatch_spline", np);
        gsFileManager::open("multipatch_spline.pvd");
    }

    return EXIT_SUCCESS;
}
