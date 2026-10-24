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
// Sample position and reference-gradient errors of Fapp versus F o Psi.
// Boundary error uses only external patch sides; these are sampled maxima.
void ComputesErrorGeometry(const gsMultiPatch<> &FF,
                           const gsMultiPatch<> &MAmapping,
                           const gsMultiPatch<> &Apmapping,
                           const int numElData,
                           real_t &maxDist,
                           real_t &Binf,
                           real_t &Ginf)
{
    const index_t dim = MAmapping.dim();
    GISMO_ENSURE(numElData > 0 && (dim == 2 || dim == 3),
                 "Use a positive sampling factor and a 2D or 3D parameter domain.");
    GISMO_ENSURE(FF.nPatches() > 0 && FF.nPatches() == MAmapping.nPatches()
                 && FF.nPatches() == Apmapping.nPatches()
                 && FF.dim() == dim && Apmapping.dim() == dim
                 && MAmapping.geoDim() == dim
                 && FF.geoDim() == Apmapping.geoDim(), "Composition dimensions must match.");
    gsStopwatch timer;
    timer.restart();

    // Keep the previous sampling resolution, independently of the fitted basis.
    gsMultiPatch<> samplingDomain;
    if (dim == 2)
        samplingDomain.addPatch(gsNurbsCreator<>::BSplineRectangle(0., 0., 1., 1.));
    else
        samplingDomain.addPatch(gsNurbsCreator<>::BSplineCube(1., 0., 0., 0.));
    gsMultiBasis<> sampling(samplingDomain, true);
    while (sampling.basis(0).size() <= samplingDomain.basis(0).size() * numElData)
        sampling.uniformRefine();
    const gsMatrix<> unitGrid = sampling.basis(0).anchors();

    maxDist = Binf = Ginf = 0.;
    for (size_t p = 0; p < FF.nPatches(); ++p)
    {
        const gsMatrix<> support = MAmapping.basis(p).support();
        gsMatrix<> grid = unitGrid;
        for (index_t d = 0; d < dim; ++d)
            grid.row(d) = support(d,0) + (support(d,1)-support(d,0)) * unitGrid.row(d).array();

        gsMatrix<> mapped = MAmapping.patch(p).eval(grid);
        const gsMatrix<> bounds = FF.basis(p).support();
        // Correct roundoff at endpoints, but do not hide a map leaving its domain.
        for (index_t d = 0; d < dim; ++d)
        {
            GISMO_ENSURE(mapped.row(d).minCoeff() >= bounds(d,0)-1e-10
                         && mapped.row(d).maxCoeff() <= bounds(d,1)+1e-10,
                         "Moving mapping leaves the initial patch parameter domain.");
            mapped.row(d) = mapped.row(d).cwiseMax(bounds(d,0)).cwiseMin(bounds(d,1));
        }
        const gsMatrix<> exact = FF.patch(p).eval(mapped);
        const gsMatrix<> approx = Apmapping.patch(p).eval(grid);
        const gsMatrix<> DF = FF.patch(p).jacobian(mapped);
        const gsMatrix<> DPsi = MAmapping.patch(p).jacobian(grid);
        const gsMatrix<> Dapprox = Apmapping.patch(p).jacobian(grid);
        GISMO_ENSURE(exact.allFinite() && approx.allFinite() && DF.allFinite()
                     && DPsi.allFinite() && Dapprox.allFinite(), "Non-finite composition values or derivatives.");
        const gsVector<> positionError = (approx-exact).colwise().norm().transpose();
        maxDist = std::max(maxDist, positionError.maxCoeff());
        for (index_t q = 0; q < grid.cols(); ++q)
            Ginf = std::max(Ginf, (Dapprox.middleCols(dim*q,dim)
                         - DF.middleCols(dim*q,dim) * DPsi.middleCols(dim*q,dim)).norm());

        // Compare corresponding boundary points, excluding internal interfaces.
        for (auto bit = FF.bBegin(); bit != FF.bEnd(); ++bit)
            if (static_cast<size_t>(bit->patch) == p)
                for (index_t q = 0; q < grid.cols(); ++q)
                    if (unitGrid(bit->direction(),q) == (bit->parameter() ? 1. : 0.))
                        Binf = std::max(Binf, positionError(q));
    }
    gsInfo << "<Error> Position_Linf = " << maxDist << ", Gradient_Linf = " << Ginf
           << ", Boundary_Linf = " << Binf << ", samples/patch = " << unitGrid.cols()
           << ", CPU time = " << timer.stop() << " s\n";
}

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot              = false;
    bool noOpen            = false;
    index_t numRefine      = 0;
    index_t numLRefine     = 1;
    index_t numRefineMAE   = 4;
    index_t numReduceMAE   = 0;
    index_t numElevateMAE  = 0;
    index_t numElevate     = 0;
    index_t maxIter        = 100;
    index_t elevDegree     = 0; // degree elevation for the composition of geometry maps
    index_t id_mp          = 4; // id of the geometry in the xml file
    index_t id_rho         = 1942; // id of the density function in the xml file
    real_t IntensityMAE    = 9.;
    real_t quadValue       = 2.0;
    bool bs_nrbs           = false;
    bool last              = true;
    bool exactGeo          = false;
    bool H1            = false;
    real_t gradientWeight = 1.;

    // gsStopwatch timer;
    // timer.restart();

    // Specify the file path
    // std::string fn("pde/Bspline_ball.xml");
    // std::string fn("volumes/GshapedVolume.xml");
    // Specify the file path
    // std::string fn("pde/circle.xml");
    // std::string fn("surfaces/egg.xml");
    // std::string fn("domain2d/lake.xml");
    // std::string fn("surfaces/cylinder.xml");
    //... multipatch case
    std::string fn("pde/annulus2d_bvp.xml");
    // std::string fn("pde/surfacepoisson_sphere_bvp.xml");
    // std::string fn("pde/surfacepoisson_torus_bvp.xml");

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
    cmd.addInt("a", "id_mp", "ID of the geometry in the XML file",
                id_mp);
    cmd.addInt("b", "id_rho", "ID of the density function in the XML file",
                id_rho);
    cmd.addInt("quRule", "Quadrature rule [1:GaussLegendre,2:GaussLobatto,3:PatchRule]",
                1);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution",
                plot);
    cmd.addSwitch("no-open", "Write visualization files without opening an external viewer", noOpen);
    cmd.addSwitch("nrb", "basis use nrubs or bspline format",
                bs_nrbs);
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement",
                last);
    cmd.addSwitch("exactGeo", "Keep the exact boundary (adapt only interfaces) and use the MMPDE solver",
                exactGeo);
    cmd.addSwitch("H1", "Use value and gradient projection for the composition", H1);
    cmd.addReal("w", "gradientWeight", "Weight of the reference-gradient error in H1 projection", gradientWeight);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Load the file
    gsFileData<> fd(fn);
    // ...
    gsMultiPatch<> mpLeft, mpPsi;// Initial geometry and the resluted adaptive mapping
    fd.getId(id_mp,mpLeft);
    gsInfo << "Loaded file " << fd.lastPath() << "\n";
    // Create a gsMultipatch and add the loaded geometry
    // gsMultiPatch<> mpLeft; mpLeft.addPatch( gsNurbsCreator<>::BSplineCube(1,0,0,0) );
    //gsMultiPatch<> mpLeft; mpLeft.addPatch( gsNurbsCreator<>::NurbsSphere(1.,0.,0.,0.));
    // mpLeft = gsNurbsCreator<>::BSplineSquareGrid(1,1,1, 0.0, 0.0);
    // mpLeft = gsNurbsCreator<>::BSplineCubeGrid(1,1,1,1.,-0.5,-0.5,-0.5);

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    mpLeft.computeTopology();

    // Right-hand side function : Analytical density function rho_1 ----- //
    // Load the file
    std::string fnDensity("pde/circle.xml");
    gsFileData<> fDensity(fnDensity);
    gsFunctionExpr<> f;
    fDensity.getId(id_rho, f);
    gsInfo<<"Density function "<< f << "\n";
    typedef gsExprAssembler<>::geometryMap geometryMap;

    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###                                  Step r* : Computes the density function
    ###                                     and the multipatch adaptive mapping from a given mesh
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    gsAdaptiveMultiPatchBuilder MAE = gsAdaptiveMultiPatchBuilder(mpLeft, numRefineMAE, maxIter, IntensityMAE, numReduceMAE, numElevateMAE, exactGeo);
    // Work with F o A on [0,1]^d, including direct composition and gradients.
    // The input XML and its original knot domains are preserved by the builder.
    mpLeft = MAE.normalizedGeometry();
    auto density        = MAE.buildAnalyticDensity(f); // build the density function (we avoid composing rho o F o Psi here)
    MAE.buildMultiPatch(density);// build the adaptive mapping

    // Show all patch densities as separate fields on one reference patch.
    gsMultiPatch<> densityDomain;
    const gsMatrix<> densityGrid = density.basis(0).anchors();
    densityDomain.addPatch(density.basis(0).interpolateData(densityGrid, densityGrid));
    gsParaviewCollection<> densityPlot("ParaviewOutput/density");
    densityPlot.newTimeStep(densityDomain);
    for (size_t p = 0; p < density.nPatches(); ++p)
        densityPlot.addField(gsField<>(densityDomain, density.patch(p), true),
                             "density_patch" + std::to_string(p + 1));
    densityPlot.saveTimeStep();
    densityPlot.save();

    // Export each adaptive reference-square mapping with its element mesh.
    for (size_t p = 0; p < MAE.MAmapping.nPatches(); ++p)
        gsWriteParaview(MAE.MAmapping.patch(p),
                        "ParaviewOutput/MAE_mapping_patch" + std::to_string(p + 1), 1000, true);

    // Keep the exact composition alive for every integration and export level.
    gsMultiPatch<> directComposition;
    for (size_t p = 0; p < mpLeft.nPatches(); ++p)
        directComposition.addPatch(gsComposedGeometry<>(MAE.MAmapping.patch(p), mpLeft.patch(p)));
    directComposition.interfaces() = mpLeft.interfaces();
    directComposition.boundaries() = mpLeft.boundaries();
    const bool embedded = mpLeft.dim() < mpLeft.geoDim();
    //...

    // //------------------------------------
    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    A.options().setReal("quA", quadValue);
    // A.options().setInt("quB", 1);
    A.options().setSwitch("SameElement",false);
    gsInfo<<"Active options:\n"<< A.options() <<"\n";

    // Elements used for numerical integration
    A.setIntegrationElements(MAE.mapping_basis);
    gsExprEvaluator<> ev(A);
    // //------------------------------------
    geometryMap G       = A.getMap(mpLeft);
    geometryMap direct  = A.getMap(directComposition);

    // h-refine each basis
    index_t numLevels = 0;
    if (last)
    {
        for (int r =0; r < numRefine; ++r){
            MAE.mapping_basis.uniformRefine();
            numLevels += 1;
        }
        numRefine = 0;
    }
    // ...
    gsVector<>  CHdferror(numRefine+1)
                , IBdrerr(numRefine+1), IPosError(numRefine+1), IGradError(numRefine+1), IL2Jerror(numRefine+1);
    gsVector<int>  DoFPDE(numRefine+1);
    for (int r=0; r<= numRefine; ++r)
    {
    if (r>0)
        MAE.mapping_basis.uniformRefine();
    numLevels += 1;
    //... some infos on the computational domain
    gsInfo << r <<"th iter:{ numElement " << MAE.mapping_basis.basis(0).numElements() << " degree " << MAE.mapping_basis.degree()
            <<" dim " <<MAE.mapping_basis.dim()<<" Geodim " << mpLeft.geoDim()
            <<"}------------------------------------------------------\n";

    DoFPDE[r]               = MAE.mapping_basis.size();
    // Embedded surfaces have a rectangular Jacobian: use sqrt(det(J^T J)).
    CHdferror[r] = embedded
        ? math::abs(ev.integral(meas(G) - meas(direct)))
        : math::abs(ev.integral(jac(G).det() - jac(direct).det()));

    //----------------------------------------------------------
    //...Interpolation of the mapping by collocation method !
    //----------------------------------------------------------
    mpPsi = H1 ? MAE.buildH1CompMultiPatch(gradientWeight, quadValue)
           : MAE.buildColCompMultiPatch();
    geometryMap PGI         = A.getMap(mpPsi);
    // ... Error using Interpolation method
    real_t maxDist = 0.;
    real_t Binf = 0.;
    real_t Ginf = 0.;
    ComputesErrorGeometry(mpLeft, MAE.MAmapping, mpPsi, 200, maxDist, Binf, Ginf);
    IL2Jerror[r] = embedded
        ? math::abs(ev.integral(meas(PGI) - meas(direct)))
        : math::abs(ev.integral(jac(PGI).det() - jac(direct).det()));
    IPosError[r]           = maxDist;
    IGradError[r]          = Ginf;
    IBdrerr[r]             = Binf;

    // Evaluate F(Psi) directly, without projecting it into a spline space.
    if (mpPsi.nPatches() != 0)
    {
        geometryMap approximation = A.getMap(mpPsi);
        for (const auto * geometry : {&directComposition, &mpPsi})
        {
            const std::string name = geometry == &directComposition ? "composition_direct" : "composition_approx";
            gsParaviewCollection<> comparison("ParaviewOutput/" + name, ev);
            comparison.options().setSwitch("elements", true);
            comparison.options().setInt("elementResolution", 8);
            comparison.options().setInt("numPoints", 10000);
            comparison.options().setInt("precision", 10);
            comparison.newTimeStep(*geometry);
            if (embedded)
            {
                comparison.addField(meas(direct), "direct_Jacobian");
                comparison.addField(meas(approximation), "approximated_Jacobian");
            }
            else
            {
                comparison.addField(jac(direct).det(), "direct_Jacobian");
                comparison.addField(jac(approximation).det(), "approximated_Jacobian");
            }
            comparison.addField((approximation - direct).norm(), "composition_error");
            comparison.addField((jac(approximation) - jac(direct)).norm(), "gradient_error");
            comparison.saveTimeStep();
            comparison.save();
        }
    }

    } // All refinement levels are complete before writing the error vectors.

    std::ofstream outFile("../build/errorGeometry_analysis.txt", std::ios::app); // Open file in append mode
    if (outFile.is_open())
    {
        outFile << "#DoF_PDE: q"<< quadValue<<"pPr"<< MAE.mapping_basis.basis(0).maxDegree()<<"pPsi"<< mpLeft.basis(0).maxDegree()-numReduceMAE <<"\n"
                << std::scientific << DoFPDE.transpose() << "\n";

        outFile << (H1 ? "#H1 projection error analysis: \n" : "#Collocation projection error analysis: \n");
        outFile << "#Position_Linf_error: \n" << std::scientific << std::setprecision(3) << IPosError.transpose() << "\n";
        outFile << "#Gradient_Linf_error: \n" << std::scientific << std::setprecision(3) << IGradError.transpose() << "\n";
        outFile << "#Boundary_Linf_error: \n" << std::scientific << std::setprecision(3) << IBdrerr.transpose() << "\n";
        outFile << "#Jacobian_integral_error: \n" << std::scientific << std::setprecision(3) << IL2Jerror.transpose() << "\n";

        // outFile << "#C_error:  "<< quadValue << ": "<< std::scientific << std::setprecision(3) << CHdferror.transpose() << "\n";
        outFile << "#-------------------------------------------------------------------------------\n"; // Optional separator for readability
        outFile.close(); // Close the file after writing
        gsInfo << "Error analysis results appended to errorGeometry_analysis.txt.\n";
    }
    else
    {
        gsInfo << "Error: Unable to open file for writing : errorGeometry_analysis.txt.\n";
    }
    //--------------------------------------------------------------------------------------------------
    //! [Error and convergence rates]
    //--------------------------------------------------------------------------------------------------
    // --- Print header ---
    gsInfo << std::setw(12) << "DoFs" << " & "
         << std::setw(13) << "Boundary"  << " & "
         << std::setw(6)  << "EOcB"      << " & "
         << std::setw(13) << "PositionInf" << " & "
         << std::setw(6)  << "EOcP"      << " & "
         << std::setw(13) << "GradientInf" << " & "
         << std::setw(6)  << "EOcG"      << " & "
         << std::setw(13) << "JacInt" << " & "
         << std::setw(6)  << "EOcJ"      << "\n";
    // --- Print table row by row ---
    gsInfo << std::string(50, '-') << (H1 ? "H1Proj\n" : "Colloc\n");
        auto orderofConvBndr = ( IBdrerr.head(numRefine).array() /
                  IBdrerr.tail(numRefine).array() ).log().transpose() / std::log(2.0);
        auto orderofConv = ( IPosError.head(numRefine).array() /
                  IPosError.tail(numRefine).array() ).log().transpose() / std::log(2.0);
        auto orderofConvGrad = ( IGradError.head(numRefine).array() /
                  IGradError.tail(numRefine).array() ).log().transpose() / std::log(2.0);
        auto orderofConvL2 = ( IL2Jerror.head(numRefine).array() /
                  IL2Jerror.tail(numRefine).array() ).log().transpose() / std::log(2.0);
    gsInfo << std::setw(12) << DoFPDE[0] << " & "
            << std::setw(12) <<std::setprecision(3)<<std::scientific<< IBdrerr[0] << "&"
            << std::setw(12) << "--" << " & "
            << std::setw(12) <<std::setprecision(3)<<std::scientific<< IPosError[0] << " & "
            << std::setw(12) << "--" << " & "
            << std::setw(12) <<std::setprecision(3)<<std::scientific<< IGradError[0] << " & "
            << std::setw(12) << "--" << " & "
            << std::setw(12) <<std::setprecision(3)<<std::scientific<< IL2Jerror[0] << "&"
            << std::setw(12) << "--" << "\\\\ \n";
    for (int i = 1; i <= numRefine; i++) {
        gsInfo << std::setw(12) << DoFPDE[i] << " & "
             << std::setw(12) <<std::setprecision(3)<<std::scientific<< IBdrerr[i] << "&"
             << std::setw(12) <<std::fixed<<std::setprecision(2)<< orderofConvBndr[i-1] << " & "
             << std::setw(12) <<std::setprecision(3)<<std::scientific<< IPosError[i] << " & "
             << std::setw(12) <<std::fixed<<std::setprecision(2)<< orderofConv[i-1] << " & "
             << std::setw(12) <<std::setprecision(3)<<std::scientific<< IGradError[i] << " & "
             << std::setw(12) <<std::fixed<<std::setprecision(2)<< orderofConvGrad[i-1] << " & "
             << std::setw(12) <<std::setprecision(3)<<std::scientific<< IL2Jerror[i] << " & "
             << std::setw(12) <<std::fixed<<std::setprecision(2)<< orderofConvL2[i-1] << "\\\\ \n";
            }
    //! [Export visualization in ParaView]
    if (plot)
    {
        gsMultiPatch<> Psi;
        // check first if the mapping is already thb-spline
        if (dynamic_cast<const gsHTensorBasis<2> *>(&mpPsi.basis(0).source()))
            Psi = mpPsi;
        else if (dynamic_cast<const gsHTensorBasis<3> *>(&mpPsi.basis(0).source()))
            Psi = mpPsi;
        else{
            gsInfo<<"The basis is not a THB-spline basis, we will use the original mapping \n";
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
        real_t adaptRefParam = 0.7;

        for (int r=0; r<numLRefine; ++r)
        {
        //! [beginRefLoop]
            gsInfo << "====== Loop " << r << " of "
                    <<numLRefine<< " ====adapt Parameter ="<< adaptRefParam << " ======" << "\n";
            // --------------- error estimation/computation ---------------
            // Get the element-wise norms.
            // Integrate density in the patch-by-patch order used for marking.
            ev.integralElWise( ff_TG.val() );
            const std::vector<real_t> eltErrs = ev.elementwise();
            GISMO_ENSURE(eltErrs.size() == basis.totalElements(),
                         "Density indicators must cover every integration element.");
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
        gsParaviewCollection<> collection("ParaviewOutput/solution", ev);
        collection.options().setSwitch("elements", true);
        collection.options().setSwitch("base64", false);
        collection.options().setInt("elementResolution", 16);
        collection.options().setInt("numPoints", 10000);
        collection.newTimeStep(Psi);
        collection.addField(ff_TG, "density function");
        if (embedded)
            collection.addField(meas(PPsi), "Jacobian function");
        else
            collection.addField(jac(PPsi).det(), "Jacobian function");
        collection.saveTimeStep();
        collection.save();

        if (!noOpen)
            gsFileManager::open("ParaviewOutput/solution.pvd");
        // gsWrite(Psi, "Psi_mapping");
        // gsInfo << "Result written in Psi_mapping.xml \n";
    }
    else
        gsInfo << "Done. Density, MAE meshes and composition comparison exported to ParaviewOutput/; re-run with --plot "
                  "to also export the solution.\n";

    return EXIT_SUCCESS;

}// end main
