/** @file r_refinement_example.cpp

    @brief Tutorial on how to use compute the Monge-Ampere mapping for adaptive mesh in two and three dimensions.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris & M. BAHARI
*/

//! [Include namespace]
#include <gismo.h>
#include <fstream>  // For file operations
#include <gsAssembler/gsAdaptiveMultiPatchBuilder.h>  // Include the new class of r_refinement

using namespace gismo;
//! [Include namespace]
// computes the projection of a composition and return a MultiPatch object :: fitting
void ComputesErrorGeometry(const gsMultiPatch<> &FF,
                           const gsMultiPatch<> &MAmapping,
                           const gsMultiPatch<> &Apmapping,
                           const int numElData,
                           double &maxDist,
                           double &Binf)
{
    gsInfo << "<Error> Compute geometric error: ";
    assert(MAmapping.dim() == 2 && "Only single-patch 2D fitting is implemented so far.");

    gsStopwatch timer;
    timer.restart();

    // --- Generate a refined sampling grid ------------------------------------
    auto corners         = MAmapping.basis(0).support();
    gsMultiPatch<> mptmp;
    mptmp.addPatch(gsNurbsCreator<>::BSplineRectangle(corners.at(0), corners.at(1), corners.at(2), corners.at(3)));
    gsMultiBasis<> Tbasis(mptmp, true);
    while (Tbasis.basis(0).size() <= mptmp.basis(0).size() * numElData)
        Tbasis.uniformRefine();
    const gsMatrix<> intGrid = Tbasis.basis(0).anchors();
    gsInfo << ": gridsize from= " << Tbasis.basis(0).size()<< " ./. ";

    // --- Evaluate mappings ---------------------------------------------------
    gsMatrix<> intVals = MAmapping.patch(0).eval(intGrid);
    gsMatrix<> JVals = MAmapping.patch(0).jacobian(intGrid);// Jacobain of square to square mapping
    intVals = intVals.cwiseMax(0).cwiseMin(1);

    const gsMatrix<> XF = FF.patch(0).eval(intVals);     // reference geometry
    const gsMatrix<> JF = FF.patch(0).jacobian(intVals);     // Jacobian of reference geometry in adapted grids by  MA mapping
    const gsMatrix<> XG = Apmapping.patch(0).eval(intGrid); // approximate geometry
    const gsMatrix<> JG = Apmapping.patch(0).jacobian(intGrid); // Jacobain of approximate geometry
    const index_t Nf = XF.cols(), Ng = XG.cols();

    // --- Initialize outputs --------------------------------------------------
    maxDist = 0.0;   // Hausdorff distance
    Binf    = 0.0;   // boundary max distance
    index_t ngrids = sqrt(intGrid.cols()); // number of grid points in one direction, assuming a square grid.
    const index_t lookAround = 50; // look around 20 points to find the closest point on the other geometry for boundary error estimation

    // --- Compute symmetric Hausdorff distance -------------------------------
    #pragma omp parallel for reduction(max:maxDist)
    for (index_t ix = 0; ix < ngrids; ++ix)
    {
    for (index_t jx = 0; jx < ngrids; ++jx)
    {
        index_t i  = ix * ngrids + jx;
        real_t minDist = std::numeric_limits<real_t>::max();
        for (index_t iy = std::max(ix-lookAround, 0); iy < std::min(ix+lookAround,ngrids); ++iy)
        {
        for (index_t jy = std::max(jx-lookAround, 0); jy < std::min(jx+lookAround,ngrids); ++jy)
        {
            index_t j  = iy * ngrids + jy;
            const real_t d = (XF.col(i) - XG.col(j)).norm();
            if (d < minDist) minDist = d;
        }
        }
        if (minDist > maxDist) maxDist = minDist;
    }
    }

    #pragma omp parallel for reduction(max:maxDist)
    for (index_t j = 0; j < Ng; ++j)
    {
        real_t minDist = std::numeric_limits<real_t>::max();
        for (index_t i = 0; i < Nf; ++i)
        {
            const real_t d = (XG.col(j) - XF.col(i)).norm();
            if (d < minDist) minDist = d;
        }
        if (minDist > maxDist) maxDist = minDist;
    }

    // --- Compute boundary only error ----------------------------------------
    #pragma omp parallel for reduction(max:Binf)
    for (index_t i = 0; i < ngrids; ++i)
    {
        // y = 0
        index_t k  = i * ngrids + 0;
        real_t minDist = std::numeric_limits<real_t>::max();
        for (index_t ii = std::max(i-lookAround, 0); ii < std::min(i+lookAround,ngrids); ++ii)
        {
        for (index_t ij = 0; ij < std::min(lookAround,ngrids); ++ij)
        {
            index_t kk  = ii * ngrids + ij;
            const real_t d = (XG.col(k) - XF.col(kk)).norm();
            if (d < minDist) minDist = d;
        }
        }
        if (minDist > Binf) Binf = minDist;
        // y = 1
        k  = i * ngrids + ngrids - 1;
        minDist = std::numeric_limits<real_t>::max();
        for (index_t ii = std::max(i-lookAround, 0); ii < std::min(i+lookAround,ngrids); ++ii)
        {
        for (index_t ij = std::max(ngrids-1-lookAround, 0); ij < ngrids; ++ij)
        {
            index_t kk  = ii * ngrids + ij;
            const real_t d = (XG.col(k) - XF.col(kk)).norm();
            if (d < minDist) minDist = d;
        }
        }
        if (minDist > Binf) Binf = minDist;
        // x = 0
        k  = 0 * ngrids + i;
        minDist = std::numeric_limits<real_t>::max();
        for (index_t ii = 0; ii < std::min(lookAround,ngrids); ++ii)
        {
        for (index_t ij = std::max(i-lookAround, 0); ij < std::min(i+lookAround,ngrids); ++ij)
        {
            index_t kk  = ii * ngrids + ij;
            const real_t d = (XG.col(k) - XF.col(kk)).norm();
            if (d < minDist) minDist = d;
        }
        }
        if (minDist > Binf) Binf = minDist;
        // x = 1
        k  = (ngrids - 1) * ngrids + i;
        minDist = std::numeric_limits<real_t>::max();
        for (index_t ii = std::max(ngrids-1-lookAround, 0); ii < ngrids; ++ii)
        {
        for (index_t ij = std::max(i-lookAround, 0); ij < std::min(i+lookAround,ngrids); ++ij)
        {
            index_t kk  = ii * ngrids + ij;
            const real_t d = (XG.col(k) - XF.col(kk)).norm();
            if (d < minDist) minDist = d;
        }
        }
        if (minDist > Binf) Binf = minDist;
    }

    double cpu = timer.stop();
    gsInfo << "Hausdorff = " << maxDist
           << ", Boundary = " << Binf
           << ", CPU time = " << cpu << " s\n";
};

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot              = false;
    index_t numRefine      = 2;
    index_t numLRefine     = 1;
    index_t numRefineMAE   = 3;
    index_t numReduceMAE   = 0;
    index_t numElevateMAE  = 0;
    index_t numElevate     = 0;
    index_t maxIter        = 50;
    index_t elevDegree     = 0; // degree elevation for the composition of geometry maps
    double IntensityMAE    = 9.;
    double quadValue       = 2.0;
    bool bs_nrbs           = false;
    bool last              = true;
    bool colloc            = false;
    bool fit               = false;
    bool L2            = true;

    // gsStopwatch timer;
    // timer.restart();

    // Specify the file path
    // std::string fn("pde/example3D.xml");
    // std::string fn("volumes/GshapedVolume.xml");
    // Specify the file path
    std::string fn("pde/bsimple.xml");
    // std::string fn("pde/infinit_plate.xml");
    // std::string fn("pde/circle.xml");
    // std::string fn("surfaces/egg.xml"); 
    // std::string fn("domain2d/lake.xml");
    // std::string fn("surfaces/cylinder.xml");

    gsCmdLine cmd("Tutorial on solving a non-linear Monge-Ampere problem.");
    cmd.addInt("i", "iter", "Maximum number of iterations for the iterative Picard", 
                maxIter);
    cmd.addInt( "e", "degreeElevation", "Number of degree Elevation steps to perform before solving (0: equalize degree in all directions)", 
                numElevate );
    cmd.addInt( "v", "elevDegree", "Number of degree Elevation steps to perform fro the composition (0: equalize degree in all directions)", 
                elevDegree );
    cmd.addInt( "r", "degreeRedution", "Number of degree Reduction steps to perform before solving MAE", 
                numReduceMAE );
    cmd.addInt( "s", "degreeMAEElevation", "Number of degree Elevation steps to perform before solving MAE", 
                numElevateMAE );
    cmd.addInt( "u", "uniformRefine", "Number of Uniform h-refinement loops",  
                numRefine );
    cmd.addInt( "m", "uniformRefineMAE", "Number of Uniform h-refinement loops for MAE mapping",  
                numRefineMAE );
    cmd.addInt( "l", "numLRefine", "Number of local h-refinement loops",  
                numLRefine );
    cmd.addReal( "f", "IntensityMAE", "Intensity of density function",  
                IntensityMAE);
    cmd.addReal("q","quadValue", "Quadrature rule number of  points in each direction", 
                quadValue);
    cmd.addString( "d", "file", "Input XML file data", 
                fn );
    cmd.addInt("quRule", "Quadrature rule [1:GaussLegendre,2:GaussLobatto,3:PatchRule]",
                1);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", 
                plot);
    cmd.addSwitch("nrb", "basis use nrubs or bspline format",
                bs_nrbs);
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement",
                last);
    cmd.addSwitch("colloc", "Compute the the compodition using collocation method",
                colloc);                
    cmd.addSwitch("fit", "Use fitting to compute the composition",
                fit);
    cmd.addSwitch("L2", "Use L2-projection to compute the composition",
                L2);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if(fit)
      bs_nrbs = true;// fitting doesn't support nurbs format
    // Load the file
    gsFileData<> fd(fn);
    // ...
    gsMultiPatch<> mpLeft, mpPsi;// Initial geometry and the resluted adaptive mapping
    fd.getId(1,mpLeft);
    gsInfo << "Loaded file " << fd.lastPath() << "\n";
    // Create a gsMultipatch and add the loaded geometry
    // gsMultiPatch<> mpLeft; mpLeft.addPatch( gsNurbsCreator<>::BSplineCube(1,0,0,0) );
    //gsMultiPatch<> mpLeft; mpLeft.addPatch( gsNurbsCreator<>::NurbsSphere(1.,0.,0.,0.));
    // mpLeft = gsNurbsCreator<>::BSplineSquareGrid(1,1,1, 0.0, 0.0);
    // mpLeft = gsNurbsCreator<>::BSplineCubeGrid(1,1,1,1.,-0.5,-0.5,-0.5);

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    mpLeft.computeTopology();
    
    // Right-hand side function : Analytical density function rho_1
    // Load the file
    gsFunctionExpr<> f;
    fd.getId(2003, f);
    gsInfo<<"Density function "<< f << "\n";

    //! [Refinement]
    gsMultiBasis<double> dbasis(mpLeft, bs_nrbs);//true: poly-splines (not NURBS)
    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis.setDegree( dbasis.maxCwiseDegree() + numElevate);
    //dbasis.degreeIncrease(numElevate);
    
    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    A.options().setReal("quA", quadValue);
    // A.options().setInt("quB", 1);
    A.options().setSwitch("SameElement",false);
    gsInfo<<"Active options:\n"<< A.options() <<"\n";

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);

    typedef gsExprAssembler<>::geometryMap geometryMap;

    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###                                  Step r* : Computes the density function
    ###                                     and the multipatch adaptive mapping from a given mesh
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    gsAdaptiveMultiPatchBuilder MAE = gsAdaptiveMultiPatchBuilder(mpLeft, numRefineMAE, maxIter, IntensityMAE, numReduceMAE, numElevateMAE);
    auto density        = MAE.buildAnalyticDensity(f); // build the density function (we avoid composing rho o F o Psi here)
    MAE.buildMultiPatch(density, 1e-8);// build the adaptive mapping
    // //------------------------------------
    geometryMap G       = A.getMap(mpLeft);
    geometryMap PP      = A.getMap(MAE.MAmapping);
    auto Cmp            = A.getCoeff(mpLeft, PP);

    // h-refine each basis
    index_t numLevels = 0;
    if (last)
    {
        for (int r =0; r < numRefine; ++r){
            dbasis.uniformRefine();
            numLevels += 1;
        }
        numRefine = 0;
    }
    // ... 
    gsVector<>  Bdrerr(numRefine+1), Hdferror(numRefine+1), L2Jerror(numRefine+1)// L2 projection errors
                , CHdferror(numRefine+1)
                , IBdrerr(numRefine+1), IHdferror(numRefine+1), IL2Jerror(numRefine+1)// Interpolation errors
                , FBdrerr(numRefine+1), FHdferror(numRefine+1), FL2Jerror(numRefine+1);// Fitting errors
    gsVector<int>  DoFPDE(numRefine+1);
    for (int r=0; r<= numRefine; ++r)
    {
    dbasis.uniformRefine();
    numLevels += 1;

    //... some infos on the computational domain
    gsInfo << r <<"th iter:{ numElement " << dbasis.basis(0).numElements() << " degree " << dbasis.degree() 
            <<" dim " <<dbasis.dim()<<" Geodim " << mpLeft.geoDim() 
            <<"}------------------------------------------------------\n";

    DoFPDE[r]               = dbasis.basis(0).size();
    CHdferror[r]            = abs(ev.integral( jac(G).det() - jac(Cmp).det()*jac(PP).det() ) );

    //----------------------------------------------------------------------
    // ... computes the composition of geometry maps L^2-projection method !
    //----------------------------------------------------------------------
    if(L2)
    {
    mpPsi           = MAE.buildCompMultiPatch(dbasis, quadValue, false);
    geometryMap GPi = A.getMap(mpPsi);
    // ... Error analysis
    double maxDist = 0.;
    double Binf    = 0.;
    ComputesErrorGeometry(mpLeft, MAE.MAmapping, mpPsi, 200, maxDist, Binf);
    L2Jerror[r]             = abs(ev.integral( jac(GPi).det() - jac(Cmp).det()*jac(PP).det() ) );
    Hdferror[r]             = maxDist;// std::abs(abs(ev.integral( meas(G)  )) - abs( ev.integral(meas(GPi)) ));
    Bdrerr[r]               = Binf;//math::sqrt(ev.integral((GPi-Cmp).sqNorm()));
    }

    if (colloc){
    //---------------------------------------------------------- 
    //...Interpolation of the mapping by collocation method !
    //----------------------------------------------------------
    mpPsi                   = MAE.buildColCompMultiPatch(dbasis);
    geometryMap PGI         = A.getMap(mpPsi);
    // ... Error using Interpolation method
    double maxDist = 0.;
    double Binf = 0.;
    ComputesErrorGeometry(mpLeft, MAE.MAmapping, mpPsi, 200, maxDist, Binf);
    IL2Jerror[r]            = abs(ev.integral( jac(PGI).det() - jac(Cmp).det()*jac(PP).det() ) );
    IHdferror[r]            = maxDist;//std::abs(abs(ev.integral( meas(G)  )) - abs( ev.integral(meas(PGI)) ));
    IBdrerr[r]              = Binf;// math::sqrt(ev.integral((PGI-Cmp).sqNorm())); 
    }

    if(fit){
    //---------------------------------------------------------- 
    //...Interpolation of the mapping by fitting method !
    //----------------------------------------------------------
    gsInfo<<"Fitting the mapping ++++++" <<numLevels<<"\n";
    mpPsi                   = MAE.buildFitCompMultiPatch(dbasis, 50, 1e-7);
    geometryMap PGF         = A.getMap(mpPsi);
    // ... Error analysis
    double maxDist = 0.;
    double Binf = 0.;
    ComputesErrorGeometry(mpLeft, MAE.MAmapping, mpPsi, 200, maxDist, Binf);
    FL2Jerror[r]            = abs(ev.integral( jac(PGF).det() - jac(Cmp).det()*jac(PP).det() ) );
    FHdferror[r]            = maxDist;//std::abs(abs(ev.integral( meas(G)  )) - abs( ev.integral(meas(PGF)) ));
    FBdrerr[r]              = Binf;// math::sqrt(ev.integral((PGF-Cmp).sqNorm())); 
    }
    }

    // Assuming DoFPDE, l2err, and h1err are gsMatrix or similar types
    std::ofstream outFile("errorGeometry_analysis.txt", std::ios::app); // Open file in append mode
    if (outFile.is_open())
    {
        outFile << "#DoF_PDE: q"<< quadValue<<"pPr"<< dbasis.basis(0).maxDegree()<<"pPsi"<< mpLeft.basis(0).maxDegree()-numReduceMAE <<"\n"
                << std::scientific << DoFPDE.transpose() << "\n";
        if (L2){
        outFile << "#L2 projection error analysis: \n";
        outFile << "#Hausdorff_error: \n" << std::scientific << std::setprecision(3) << Hdferror.transpose() << "\n";
        outFile << "#Boundary_error: \n" << std::scientific << std::setprecision(3) << Bdrerr.transpose() << "\n";
        outFile << "#L2_error: \n" << std::scientific << std::setprecision(3) << L2Jerror.transpose() << "\n";
        }
        if (colloc){
        outFile << "#Collocation projection error analysis: \n";
        outFile << "#Hausdorff_error: \n" << std::scientific << std::setprecision(3) << IHdferror.transpose() << "\n";
        outFile << "#Boundary_error: \n" << std::scientific << std::setprecision(3) << IBdrerr.transpose() << "\n";
        outFile << "#L2_error: \n" << std::scientific << std::setprecision(3) << IL2Jerror.transpose() << "\n";
        }
        if (fit){
        outFile << "#Fitting projection error analysis: \n";
        outFile << "#Hausdorff_error: \n" << std::scientific << std::setprecision(3) << FHdferror.transpose() << "\n";
        outFile << "#Boundary_error: \n" << std::scientific << std::setprecision(3) << FBdrerr.transpose() << "\n";
        outFile << "#L2_error: \n" << std::scientific << std::setprecision(3) << FL2Jerror.transpose() << "\n";
        }
        // outFile << "#C_error:  "<< quadValue << ": "<< std::scientific << std::setprecision(3) << CHdferror.transpose() << "\n";
        outFile << "#-------------------------------------------------------------------------------\n"; // Optional separator for readability
        outFile.close(); // Close the file after writing
        gsInfo << "Error analysis results appended to errorGeometry_analysis.txt.\n";
    }
    else
    {
        gsInfo << "Error: Unable to open file for writing : error_analysis.txt.\n";
    }
    //--------------------------------------------------------------------------------------------------
    //! [Error and convergence rates]
    //--------------------------------------------------------------------------------------------------    
    // --- Print header ---
    gsInfo << std::setw(12) << "DoFs" << " & "
         << std::setw(13) << "Boundary"  << " & "
         << std::setw(6)  << "EOcB"      << " & "
         << std::setw(13) << "Hausdorff"     << " & "
         << std::setw(6)  << "EOcH"      << " & "
         << std::setw(13) << "L2" << " & "
         << std::setw(6)  << "EOcL2"      << "\n";
    // --- Print table row by row ---
    if (L2){
    gsInfo << std::string(50, '-') << "L2Proj\n";
        auto orderofConvBndr = ( Bdrerr.head(numRefine).array() /
                  Bdrerr.tail(numRefine).array() ).log().transpose() / std::log(2.0);
        auto orderofConv = ( Hdferror.head(numRefine).array() /
                  Hdferror.tail(numRefine).array() ).log().transpose() / std::log(2.0);
        auto orderofConvL2 = ( L2Jerror.head(numRefine).array() /
                  L2Jerror.tail(numRefine).array() ).log().transpose() / std::log(2.0);
    gsInfo << std::setw(12) << DoFPDE[0] << " & "
            << std::setw(12) <<std::setprecision(3)<<std::scientific<< Bdrerr[0] << "&"
            << std::setw(12) << "--" << " & "
            << std::setw(12) <<std::setprecision(3)<<std::scientific<< Hdferror[0] << " & "
            << std::setw(12) << "--" << " & "
            << std::setw(12) <<std::setprecision(3)<<std::scientific<< L2Jerror[0] << "&"
            << std::setw(12) << "--" << "\\\\ \n";
    for (int i = 1; i <= numRefine; i++) {
        gsInfo << std::setw(12) << DoFPDE[i] << " & "
             << std::setw(12) <<std::setprecision(3)<<std::scientific<< Bdrerr[i] << "&"
             << std::setw(12) <<std::fixed<<std::setprecision(2)<< orderofConvBndr[i-1] << " & "
             << std::setw(12) <<std::setprecision(3)<<std::scientific<< Hdferror[i] << " & "
             << std::setw(12) <<std::fixed<<std::setprecision(2)<< orderofConv[i-1] << " & "
            << std::setw(12) <<std::setprecision(3)<<std::scientific<< L2Jerror[i] << " & "
             << std::setw(12) <<std::fixed<<std::setprecision(2)<< orderofConvL2[i-1] << "\\\\ \n";
            }
    }
    if (colloc){
    gsInfo << std::string(50, '-') << "Colloc\n";
        auto orderofConvBndr = ( IBdrerr.head(numRefine).array() /
                  IBdrerr.tail(numRefine).array() ).log().transpose() / std::log(2.0);
        auto orderofConv = ( IHdferror.head(numRefine).array() /
                  IHdferror.tail(numRefine).array() ).log().transpose() / std::log(2.0);
        auto orderofConvL2 = ( IL2Jerror.head(numRefine).array() /
                  IL2Jerror.tail(numRefine).array() ).log().transpose() / std::log(2.0);
    gsInfo << std::setw(12) << DoFPDE[0] << " & "
            << std::setw(12) <<std::setprecision(3)<<std::scientific<< IBdrerr[0] << "&"
            << std::setw(12) << "--" << " & "
            << std::setw(12) <<std::setprecision(3)<<std::scientific<< IHdferror[0] << " & "
            << std::setw(12) << "--" << " & "
            << std::setw(12) <<std::setprecision(3)<<std::scientific<< IL2Jerror[0] << "&"
            << std::setw(12) << "--" << "\\\\ \n";
    for (int i = 1; i <= numRefine; i++) {
        gsInfo << std::setw(12) << DoFPDE[i] << " & "
             << std::setw(12) <<std::setprecision(3)<<std::scientific<< IBdrerr[i] << "&"
             << std::setw(12) <<std::fixed<<std::setprecision(2)<< orderofConvBndr[i-1] << " & "
             << std::setw(12) <<std::setprecision(3)<<std::scientific<< IHdferror[i] << " & "
             << std::setw(12) <<std::fixed<<std::setprecision(2)<< orderofConv[i-1] << " & "
             << std::setw(12) <<std::setprecision(3)<<std::scientific<< IL2Jerror[i] << " & "
             << std::setw(12) <<std::fixed<<std::setprecision(2)<< orderofConvL2[i-1] << "\\\\ \n";
            }
    }
    if (fit){
    gsInfo << std::string(50, '-') << "Fit\n";
        auto orderofConvBndr = ( FBdrerr.head(numRefine).array() /
                  FBdrerr.tail(numRefine).array() ).log().transpose() / std::log(2.0);
        auto orderofConv = ( FHdferror.head(numRefine).array() /
                  FHdferror.tail(numRefine).array() ).log().transpose() / std::log(2.0);
        auto orderofConvL2 = ( FL2Jerror.head(numRefine).array() /
                  FL2Jerror.tail(numRefine).array() ).log().transpose() / std::log(2.0);
    gsInfo << std::setw(12) << DoFPDE[0] << " & "
            << std::setw(12) <<std::setprecision(3)<<std::scientific<< FBdrerr[0] << "&"
            << std::setw(12) << "--" << " & "
            << std::setw(12) <<std::setprecision(3)<<std::scientific<< FHdferror[0] << " & "
            << std::setw(12) << "--" << " & "
            << std::setw(12) <<std::setprecision(3)<<std::scientific<< FL2Jerror[0] << "&" 
            << std::setw(12) << "--" << "\\\\ \n";
    for (int i = 1; i <= numRefine; i++) {
        gsInfo << std::setw(12) << DoFPDE[i] << " & "
             << std::setw(12) <<std::setprecision(3)<<std::scientific<< FBdrerr[i] << "&"
             << std::setw(12) <<std::fixed<<std::setprecision(2)<< orderofConvBndr[i-1] << " & "
             << std::setw(12) <<std::setprecision(3)<<std::scientific<< FHdferror[i] << " & "
             << std::setw(12) <<std::fixed<<std::setprecision(2)<< orderofConv[i-1] << " & "
             << std::setw(12) <<std::setprecision(3)<<std::scientific<< FL2Jerror[i] << "&"
             << std::setw(12) <<std::fixed<<std::setprecision(2)<< orderofConvL2[i-1] << "\\\\ \n";
            }
    }
    //! [Export visualization in ParaView] 
    if (plot)
    {
        gsMultiPatch<> Psi;
        if (fit){ // already in THB format
            Psi = mpPsi;
        }
        else{
        if ( mpPsi.basis(0).weights().any()){
            gsInfo<<"Rational mapping \n";
        if (mpLeft.dim()== 3){
        for(size_t i =0; i<mpPsi.nPatches(); ++i)
            Psi.addPatch(gsRationalTHBSpline<3>( dynamic_cast<const gsTensorNurbs<3>&>(mpPsi.patch(i)) ));
        }
        else{
        for(size_t i =0; i<mpPsi.nPatches(); ++i)
            Psi.addPatch(gsRationalTHBSpline<2>( dynamic_cast<const gsTensorNurbs<2>&>(mpPsi.patch(i)) ));
        }
        }
        else{
            gsInfo<<"nonRational mapping \n";
        if (mpLeft.dim()== 3){
        for(size_t i =0; i<mpPsi.nPatches(); ++i)
            Psi.addPatch(gsTHBSpline<3>( dynamic_cast<const gsTensorBSpline<3>&>(mpPsi.patch(i)) ));
        }
        else{
        for(size_t i =0; i<mpPsi.nPatches(); ++i)
            Psi.addPatch(gsTHBSpline<2>( dynamic_cast<const gsTensorBSpline<2>&>(mpPsi.patch(i)) ));            
        }
        }
        }
        Psi.addAutoBoundaries();
        Psi.computeTopology();

        //Psi.uniformRefine();
        gsMultiBasis<> basis(Psi, true);//true: poly-splines (not NURBS)
        // Elements used for numerical integration
        A.setIntegrationElements(basis);
        gsExprEvaluator<> ev(A);
        // ...
        geometryMap PPsi = A.getMap(Psi);
        auto ff_TG       = A.getCoeff(f, PPsi);
        // --------------- adaptive refinement ---------------
        // Specify cell-marking strategy...
        MarkingStrategy adaptRefCrit = PUCA;
        //MarkingStrategy adaptRefCrit = GARU;
        //MarkingStrategy adaptRefCrit = errorFraction;
        real_t adaptRefParam = 0.8;

        for (int r=0; r<=numLRefine; ++r)
        {
        //! [beginRefLoop]
            gsInfo << "====== Loop " << r << " of "
                    <<numLRefine<< " ====adapt Parameter ="<< adaptRefParam << " ======" << "\n";
            // --------------- error estimation/computation ---------------
            // Get the element-wise norms.
            ev.integralElWise( ff_TG.val() );
            //ev.integralElWise( 1/jac(GPi).absDet() );

            const std::vector<real_t> eltErrs  = ev.elementwise();
            //! [errorComputation]

            //! [adaptRefinementPart]
            // Mark elements for refinement, based on the computed local errors and
            // the refinement-criterion and -parameter.
            std::vector<bool> elMarked( eltErrs.size() );
            gsMarkElementsForRef( eltErrs, adaptRefCrit, adaptRefParam, elMarked);
            gsInfo <<"Marked "<< std::count(elMarked.begin(), elMarked.end(), true) <<" elements.\n";
            // Refine the marked elements with a 1-ring of cells around marked elements
            gsRefineMarkedElements( basis, elMarked, 1);
            gsRefineMarkedElements( Psi, elMarked, 1);
        }

        //::::::::::::::::::::      end       :::::::::::::::::::::::::
        gsInfo<<"Plotting in Paraview...\n";
        gsParaviewCollection collection("ParaviewOutput/solution", &ev);
        collection.options().setSwitch("plotElements", true);
        collection.options().setSwitch("base64", false);
        collection.options().setInt("plotElements.resolution", 16);
        collection.options().setInt("numPoints", 10000);
        collection.newTimeStep(&Psi);
        collection.addField(ff_TG, "density function");
        collection.addField(jac(PPsi).det(), "Jacobian function");
        collection.saveTimeStep();
        collection.save();

        gsFileManager::open("ParaviewOutput/solution.pvd");
        // gsWrite(Psi, "Psi_mapping");
        // gsInfo << "Result written in Psi_mapping.xml \n";
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";

    return EXIT_SUCCESS;


}// end main