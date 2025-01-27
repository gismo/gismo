/** @file example_shell3D_CC.cpp

    @brief Analysis of a pinned flat plate

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Authors: H.M.Verhelst, C. Chianese
*/

#include <gismo.h>
#ifdef gsKLShell_ENABLED
#include <gsKLShell/gsKLShell.h>
#endif
#include <gsAssembler/gsEmbeddingUtils.h>

using namespace gismo;

#ifdef gsKLShell_ENABLED
// Choose among various shell examples, default = Thin Plate
int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot  = false;
    bool stress= false;
    bool nonlinear = false;
    index_t numRefine  = 1;
    index_t numElevate = 1;

    // Shell material properties
    real_t E_modulus = 1.0;
    real_t PoissonRatio = 0.30;
    real_t Density = 1.0;
    real_t thickness = 1.0;

    // Beam material properties
    real_t E_modulus_b = 1.0; // [MPa]
    real_t PoissonRatio_b = 0.30;
    real_t thickness_b = 1; // [mm]
    real_t height_b = 1; // [mm]
    auto G_modulus_b = 0.5 * E_modulus_b / (1+PoissonRatio_b);
    auto area_b = height_b * thickness_b;
    auto inertiamin_b = (height_b * pow(thickness_b,3))/12; // minimum principal inertia moment of beam cross-section
    auto inertiamax_b = (thickness_b * pow(height_b,3))/12; // maximum principal inertia moment of beam cross-section
    auto inertiap_b = inertiamin_b + inertiamax_b; // polar inertia moment of beam cross-section

    gsCmdLine cmd("3D shell example.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("stress", "Create a ParaView visualization file with the stresses", stress);
    cmd.addSwitch( "nl", "Nonlinear analysis", nonlinear );

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Define material parameters and geometry]
    gsMultiPatch<> mp_surf;
    gsMultiPatch<> mp_surf_def;
    gsMultiPatch<> mp_curve;

    mp_surf.addPatch( gsNurbsCreator<>::BSplineSquare(1) );
    mp_surf.embed(3);
    mp_surf.addAutoBoundaries();

    // Make a BSpline curve within the surface source domain
    gsKnotVector<> kv_c(0, 1, 2, 3); //start,end,interior knots, start/end multiplicity
    gsMatrix<> coefs_c(5, 2); //u,v;..
    coefs_c << 0, 0,
             0.2618, 0.053,
             0.5, 0.5,
             0.738, 0.9465,
             1,1;

    gsBSpline<> curve( kv_c, give(coefs_c));
    mp_curve.addPatch(curve);
    //! [Define material parameters and geometry per example]

    //! [Refine and elevate]
    // p-refine
    if (numElevate!=0)
        mp_surf.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        mp_surf.uniformRefine();

    mp_surf_def = mp_surf;
    gsWriteParaview<>(mp_surf_def, "mp_surf", 1000, true, true);
    //! [Refine and elevate]

    gsMultiBasis<> dbasis_surf(mp_surf);
    gsMultiBasis<> dbasis_curve(mp_curve);
    gsInfo << "Patches: "<< mp_surf.nPatches() <<", degree: "<< dbasis_surf.minCwiseDegree() <<"\n";
    gsInfo << dbasis_surf.basis(0)<<"\n";

    //! [Set boundary conditions and loads]
    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp_surf);

    // Point loads
    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();
    gsVector<> point(2);
    gsVector<> load(3);
    point<< 0.5,0.5;
    load << 0,0,1;
    pLoads.addLoad(point, load, 0);

    // Pressure over the whole plate
    gsVector<> tmp(3);
    tmp << 0,0,-1;
    gsConstantFunction<> force(tmp,3);

    //gsVector<> neu(3);
    // neu << 0, 1, 0;
    // gsConstantFunction<> neuData(neu,3);
    // bc.addCondition(boundary::north, condition_type::neumann, &neuData );

    for (index_t i=0; i!=3; ++i)
    {
        bc.addCondition(0,boundary::north, condition_type::dirichlet, 0, 0, false, i);
        bc.addCondition(0,boundary::east, condition_type::dirichlet, 0, 0, false,  i);
        bc.addCondition(0,boundary::west, condition_type::dirichlet, 0, 0, false,  i);
        bc.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0, false, i);
    }
    //! [Set boundary conditions and loads]

    //! [Make material functions]
    // Linear isotropic material model
    gsFunctionExpr<> t(std::to_string(thickness),3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsFunctionExpr<> rho(std::to_string(Density),3);

    std::vector<gsFunctionSet<>*> parameters;
    parameters.resize(2);
    parameters[0] = &E;
    parameters[1] = &nu;
    //! [Make material functions]

    //! [Make assembler]
    gsMaterialMatrixBase<real_t>* materialMatrix;
    gsOptionList options;
    //options.addInt("Material","Linear Isotropic",0);
    //options.addInt("Implementation","Analytical",1);
    materialMatrix = getMaterialMatrix<3,real_t>(mp_surf,t,parameters,rho,options);
    //materialMatrix->options().setInt("TensionField",0);

    // Construct the gsThinShellAssembler
    gsThinShellAssemblerBase<real_t>* assembler;
    assembler = new gsThinShellAssembler<3, real_t,true>(mp_surf,dbasis_surf,bc,force,materialMatrix);
    //assembler->options().setInt("Continuity",0);
    assembler->setPointLoads(pLoads);
    //! [Make assembler]

    ThinShellAssemblerStatus status;
    gsInfo<<"Assembling shell system\n";
    status = assembler->assemble();
    GISMO_ENSURE(status==ThinShellAssemblerStatus::Success,"Assembly failed");

    //! [Assemble linear part]
    gsSparseMatrix<> matrix = assembler->matrix();
    gsVector<> vector = assembler->rhs();
    gsInfo<<"Shell assembly done\n";

    gsInfo <<"Assembling beam system\n";
    // Declare the expression assembler
    gsExprAssembler<> exprAssembler = assembler->assembler();
    exprAssembler.cleanUp();
    // gsExprAssembler<> exprAssembler(1,1);

    // Register expressions inside assembler
    auto G_curve = exprAssembler.getMap(mp_curve);
    auto G_surf  = exprAssembler.getMap(mp_surf);
    auto defG_surf = exprAssembler.getMap(mp_surf); //NOTE: deformed surface must be set as solution field !
    auto u_surf = exprAssembler.trialSpace(0);
    // auto u_surf = exprAssembler.getSpace(assembler->getSpaceBasis(),mp_surf.targetDim(),0);

    // Expression to assemble
    auto eps = 0.5 * (ctv(defG_surf, G_curve).adj().tr()*ctv(defG_surf, G_curve) - ctv(G_surf, G_curve).adj().tr()*ctv(G_surf, G_curve));
    auto eps_der = ctv_var1(u_surf,defG_surf,G_curve) * ctv(defG_surf, G_curve);
    auto eps_der2 = ctv_var1(u_surf,defG_surf,G_curve)*ctv_var1(u_surf,defG_surf,G_curve).tr();
    auto expr = E_modulus_b/pow(ctv(G_surf, G_curve).norm(),4) * (area_b * eps_der * eps_der.tr()  + area_b * eps.val() * eps_der2);

    // Computes Dirichlet values and eliminates Dirichlet DoFs
    u_surf.setup(bc, dirichlet::interpolation, 0); // Is this needed?

    // exprAssembler.initSystem();
    index_t N = exprAssembler.numDofs();
    gsSparseMatrix<> matrix_curve(N,N);

    // Setup domain iterator
    const gsBasis<> & basis = dbasis_curve.basis(0);
    gsBasis<>::domainIter domIt = basis.makeDomainIterator();

    // Define quadrature rule
    typename gsQuadRule<>::uPtr QuRule = gsQuadrature::getPtr(basis,assembler->options().getGroup("ExprAssembler"));

    // Loop over the curve elements
    gsMatrix<> quPointsCurve, quPointsSurface, quPointsPhysical;
    gsVector<> quWeights;

    // Make expression data


    gsExprEvaluator<> exprEvaluator(exprAssembler);

    // const gsExprHelper<real_t>::Ptr exprData = exprAssembler.exprData();
    // Parse the expression
    // exprData->parse(expr);
    for (; domIt->good(); domIt->next() )
    {
        // Map the Quadrature rule to the element
        QuRule->mapTo( domIt->lowerCorner(), domIt->upperCorner(), quPointsCurve, quWeights);

        // MAYBE NOT NEEDED
        // Map the quadrature points to the surface
        mp_curve.patch(0).eval_into(quPointsCurve, quPointsSurface);
        //mp_surf.patch(0).eval_into(quPointsSurface, quPointsPhysical);
        // gsWriteParaviewPoints(quPointsPhysical,"quPointsPhysical");

        gsDebug<<"Curve integration points: "<<quPointsCurve<<"\n";
        //gsDebug<<"Surface integration points: "<<quPointsSurface<<"\n";

        // Loop over quadrature points
        for (index_t k = 0; k != quWeights.rows(); ++k)
        {

            // Precompute the expression
            // exprData->points() = quPointsCurve.col(k);
            // exprData->points() = quPointsSurface.col(k);
            // exprData->weights() = quWeights;
            // exprData->precompute(0); // updates the actives for every element
            // Create a local element matrix
            // gsMatrix<> localMat = quWeights[k] * expr.eval(0);


            gsMatrix<> evalMat = exprEvaluator.eval(expr,quPointsCurve.col(k));
            gsMatrix<> localMat = quWeights[k] * evalMat;

            // Push the local element matrix inside the big system
            const expr::gsFeSpace<real_t> & v = expr.rowVar();
            const expr::gsFeSpace<real_t> & u = expr.colVar();
            const index_t rd                  = v.dim();//row
            const index_t cd                  = u.dim();//col
            const gsDofMapper  & rowMap       = v.mapper();
            const gsDofMapper  & colMap       = u.mapper();
            gsMatrix<index_t>   colInd0       = assembler->getSpaceBasis().basis(0).active(quPointsSurface.col(k));
            gsMatrix<index_t>   rowInd0       = assembler->getSpaceBasis().basis(0).active(quPointsSurface.col(k));
            // gsMatrix<index_t> & colInd0       = const_cast<gsMatrix<index_t>&>(u.data().actives);// NOT AVAILABLE, SINCE U_SURF IS NOT PARSED BY OUR EXPRESSIONS
            // gsMatrix<index_t> & rowInd0       = const_cast<gsMatrix<index_t>&>(v.data().actives);// NOT AVAILABLE, SINCE U_SURF IS NOT PARSED BY OUR EXPRESSIONS

            // Push
            for (index_t r = 0; r != rd; ++r)
            {
                const index_t rls = r * rowInd0.rows();     //local stride
                for (index_t i = 0; i != rowInd0.rows(); ++i)
                {
                    // const index_t ii = rowMap.index(rowInd0.at(i),v.data().patchId,r); //N_i
                    const index_t ii = rowMap.index(rowInd0.at(i),0,r); //N_i
                    if ( rowMap.is_free_index(ii) )
                    {
                        for (index_t c = 0; c != cd; ++c)
                        {
                            const index_t cls = c * colInd0.rows();     //local stride
                            for (index_t j = 0; j != colInd0.rows(); ++j)
                            {
                                if ( 0 == localMat(rls+i,cls+j) ) continue;

                                // const index_t jj = colMap.index(colInd0.at(j),u.data().patchId,c); // N_j
                                const index_t jj = colMap.index(colInd0.at(j),0,c); // N_j
                                if ( colMap.is_free_index(jj) )
                                {
                                    matrix_curve.coeffRef(ii, jj) += localMat(rls+i,cls+j);
                                }
                            }
                        }
                    }
                }
            }

        }

        // // Create a local element matrix
        // gsMatrix<> localMat;
        // // Integrate the local element matrix
        // localMat = quWeights[0] * expr.eval(0);
        // for (index_t k = 1; k != quWeights.rows(); ++k)
        //     localMat += quWeights[k] * expr.eval(k);

        // gsDebugVar(localMat);

    }
    gsInfo<<"Assembly done\n";
    //! [Assemble linear part]

    matrix += matrix_curve;


    //! [Solve linear problem]
    gsInfo<<"Solving system with "<<assembler->numDofs()<<" DoFs\n";
    gsVector<> solVector;
    gsSparseSolver<>::CGDiagonal solver;
    solver.compute( matrix );
    solVector = solver.solve(vector);
    //! [Solve linear problem]

    //! [Construct and evaluate solution]
    mp_surf_def = assembler->constructSolution(solVector);

    gsMultiPatch<> deformation = mp_surf_def;
    for (size_t k = 0; k != mp_surf_def.nPatches(); ++k)
        deformation.patch(k).coefs() -= mp_surf.patch(k).coefs();

    gsInfo <<"Maximum deformation coef: "
           << deformation.patch(0).coefs().colwise().maxCoeff() <<".\n";
    gsInfo <<"Minimum deformation coef: "
           << deformation.patch(0).coefs().colwise().minCoeff() <<".\n";
    //! [Construct and evaluate solution]

    // ! [Export visualization in ParaView]
    if (plot)
    {
        gsField<> solField(mp_surf_def, deformation);
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( solField, "Deformation", 1000, true);
        gsWriteParaview<>( mp_surf_def, "mp_surf_def", 1000, true,true);
    }

    delete assembler;
    delete materialMatrix;
    return EXIT_SUCCESS;
}
#else
int main(int argc, char *argv[])
{
    gsInfo << "This example requires gsKLShell. Please reconfigure with gsKLShell_ENABLED=ON.\n";
    return EXIT_FAILURE;
}
#endif