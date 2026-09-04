/** @file gsAdaptiveMultiPatchBuilder.cpp

    @brief Provides generic routines for adaptive refinement.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. BAHARI
*/

#include <gismo.h>  // Include necessary GISMO headers
#include "gsAdaptiveMultiPatchBuilder.h"

typedef gsExprAssembler<>::geometryMap geometryMap;
typedef gsExprAssembler<>::variable    variable;
typedef gsExprAssembler<>::space       space;
typedef gsExprAssembler<>::solution    solution;
gsSparseSolver<>::CGDiagonal solver; // ... for the composition


//+++++++++++++++++++++++++++
// Constructor implementation
gsAdaptiveMultiPatchBuilder::gsAdaptiveMultiPatchBuilder(const gsMultiPatch<> mapping,
                            index_t numRefine,
                            index_t maxIter,
                            double IntensityMAE,
                            index_t numReduce,
                            index_t numElevate)
{
    gsInfo<<"\n <>r-refinement ";
    //-------------------------------------------------------------------------------------
    // Build a (B-spline) multi-basis for the Monge–Ampère solver
    gsMultiBasis<> dbasis(mapping, true);
    //... refine basis for convergence 
    for (int r=0; r<numRefine; ++r)
        dbasis.uniformRefine();
    // Reduce degree if possible while maintaining minimum degree of 1
    if (numReduce > 0){
        int reduceDegree = std::min(dbasis.degree()-1, numReduce);
        dbasis.degreeDecrease(reduceDegree);
    }
    else if (numElevate > 0){
        // Elevate degree if possible
        dbasis.degreeElevate(numElevate);
    }
    // Store input parameters
    this->m_basis        = dbasis;
    this->m_maxIter      = maxIter;
    this->m_IntensityMAE = IntensityMAE;
    this->DoFs           = m_basis.size();
    this->initial_mapping= mapping; 

    gsInfo << "nb patches" << mapping.nPatches() << "Using B-splines of degree " << dbasis.degree() << " DoFs ";

    //--------------------------------------------------
    // Create parametric identity mapping
    auto corners         = dbasis.basis(0).support();
    gsMultiPatch<> mp;
    mp.addPatch(gsNurbsCreator<>::BSplineRectangle(corners.at(0), corners.at(1), corners.at(2), corners.at(3)));
    // BSplineSquareGrid(1,1,corners.at(2), corners.at(0), corners.at(1));
    neumann_id = gsFunctionExpr<>("x","y",2);
    if (dbasis.dim() == 3){
        mp  = gsNurbsCreator<>::BSplineCubeGrid(1,1,1,corners.at(3),corners.at(0)-0.5,corners.at(1)-0.5,corners.at(2)-0.5);
        // Manufactured identity mapping
        neumann_id = gsFunctionExpr<>("x","y","z",3);
    }
    // Compute topology for interfaces and boundaries
    mp.computeTopology();
    this->identity_mp    = mp;

    //--------------------------------------------------
    // Initialize boundary conditions for Monge-Ampere solver
    this->bc_mae.setGeoMap(identity_mp);
    for ( gsMultiPatch<>::const_biterator
            bit = identity_mp.bBegin(); bit != identity_mp.bEnd(); ++bit)
    {
        this->bc_mae.addCondition( *bit, condition_type::neumann, &neumann_id );
    }

    //::::::::::::::::::::      Poisson fast diagonalization solver         :::::::::::::::::::::::::
    gsExprAssembler<> A(1,1);
    gsPatchPreconditionersCreator<double>::Poisson_FastDiag Poisson(dbasis.basis(0), this->bc_mae, A.options(), 1e-6);  

    // Store input parameters
    this->Poisson        = Poisson;
    gsInfo<<this->DoFs <<"<> \n";
}


/** \brief Refine every basis uniformly.
 * \param bnumRefine Number of Uniform h-refinement times.
 */
void gsAdaptiveMultiPatchBuilder::uniformRefine(const index_t numRefine)
{
    //--------------------------------------------------
    //refine basis and mapping (as we update mapping coefs/the basis)
    for( index_t i=0; i< numRefine; ++i){
    this->m_basis.uniformRefine();
    this->MAmapping.uniformRefine();
    }

    //------------
    // update DoFs
    this->DoFs      = m_basis.size();
    
    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    
    //::::::::::::::::::::      Poisson fast diagonalization solver         :::::::::::::::::::::::::
    gsPatchPreconditionersCreator<double>::Poisson_FastDiag Poisson(this->m_basis.basis(0), bc_mae, A.options(), 1e-6);  
    this->Poisson   = Poisson;
    gsInfo << "<>"<< this->DoFs <<" DoFs after uniRefine <>\n";
}


/** \brief Build and return a density as a MultiPatch object from analytical function 
 * \remark (we avoid three compositions (r o F o Psi) here to be r o Psi)
 * \param f density function defined on the physical domaine.
 */
gsMultiPatch<> gsAdaptiveMultiPatchBuilder::buildAnalyticDensity(const gsFunctionExpr<>   &f) const 
{
    gsInfo<<"<>density function: ";

    //! Assemnbler
    gsExprAssembler<> A(1,1);

    // Elements used for numerical integration
    A.setIntegrationElements(this->m_basis);

    //.. Evaluator
    gsExprEvaluator<> ev(A);

    // Set the discretization space
    space u             = A.getSpace(this->m_basis);

    // Set the Target geometry map
    geometryMap GLeft   = A.getMap(this->initial_mapping);

    // Set the source term with respect to target geometry
    auto ff             = A.getCoeff(f, GLeft);

    // Solution vector and solution variable
    gsMatrix<> densityVector;

    //u.setup(bc_mae, dirichlet::l2Projection, 0);
    A.initSystem();
    A.assemble(u* ff.val()); //rhs vector

    // prjection L2 of the composition
    densityVector        = this->Poisson.L2ProjectScalar(A.rhs());
    
    // save density funcrtion as multipatch
    solution density_sol = A.getSolution(u, densityVector);
    gsMultiPatch<> density;
    density_sol.extract(density);
    gsInfo<< densityVector.maxCoeff()<< "/"<<densityVector.minCoeff()<<"<>\n";
    return  density;
}


/**
 * \brief Build and return a density as a MultiPatch object from a solution vector using local h-refinement strategies.
 *
 * \remark We construct a new basis for the density function in order to control the mesh distribution over the computational domain.
 * If a given density is defined on an adaptive mesh, it must first be projected to a uniform mesh 
 * since the initial mapping will be used in the composition.  (r o F o Psi) to (r o F)
 * \param Hbasis Basis functions where the PDE solution is defined.
 * \param elMarked Vector containing the marked elements.
 * \param setRhogrid If 0, the density is reset to zero before adding the error distribution; otherwise, the error is accumulated.
 * \param setRhoZero Density initialization flag: if 0, the density is set to zero; otherwise, it is set to a user-defined value (fixed here to 0.75).
 */
gsMultiPatch<> gsAdaptiveMultiPatchBuilder::buildDensity(const gsMultiBasis<> Hbasis, const  std::vector<bool> elMarked, const index_t setRhogrid, const index_t setRhoZero) const 
{   
    gsInfo<<"<>density function";

    // ... basis_0 of degree 0 to represent error as a piecewise constant function
    gsMultiBasis<> basis_0 (identity_mp, true); // make a copy of basis before adaptive refinement
    basis_0.uniformRefine(setRhogrid);
    while (basis_0.size() <= std::min( this->m_basis.size(), Hbasis.size())) // refine until having enough resolution for error representation
        basis_0.uniformRefine(); // refine to have enough resolution for error representation

    // ... We want each element to be reprensted by one basis for all patches
    for (size_t pn=0; pn < this->initial_mapping.nPatches(); ++pn ) 
    {
    for( index_t i_dir=0; i_dir<basis_0.dim(); ++i_dir){
       basis_0.basis(pn).degreeDecrease(basis_0.basis(pn).degree(i_dir),i_dir);
       }
    }
    gsInfo << " Using degree "<< basis_0.basis(0).size() <<"DoFs for Piecewise Cst rho";
    // ...
    gsExprAssembler<> A_0(1,1);
    // Elements used for numerical integration
    A_0.setIntegrationElements(basis_0);    
    // Set the discretization space
    space u_0          = A_0.getSpace(basis_0);
    A_0.initSystem();

    // if (this->errorVector.size() < A_0.rhs().size() || setRhoLevel == 0){
    //     gsInfo << "~buildfromzero~";
    //     this->errorVector.resize(A_0.rhs().size());
    //     // ... initialize to zero
    //     this->errorVector.setZero();
    // }
    gsInfo << "~";
    this->errorVector.resize(A_0.rhs().size());
    gismo::gsMatrix<> ElmAccuNmb;
    ElmAccuNmb.resize(A_0.rhs().size());
    // ... initialize to zero
    this->errorVector.setZero();
    ElmAccuNmb.setOnes();
    solution error_sol = A_0.getSolution(u_0, errorVector);

    //------------------------------------------------------
    // piecewise density construction from error distribution 
    // globalCount: counter for the current global element index
    int globalCount = 0;
    #pragma omp parallel for
    for (size_t pn=0; pn < this->initial_mapping.nPatches(); ++pn )// for all patches
    {
        // for all elements in patch pn
        typename gsBasis<>::domainIter domIt =  // add patchInd to domainiter ?
            Hbasis.basis(pn).domain()->beginAll();
        typename gsBasis<>::domainIter domItEnd =  // add patchInd to domainiter ?
            Hbasis.basis(pn).domain()->endAll();
        #pragma omp parallel for
        for (; domIt<domItEnd; ++domIt )
        {
            if( elMarked[ globalCount++ ] ){ // refine this element ?
                // element index in the basis_0
                auto gIndex = basis_0.basis(pn).elementIndex(domIt.centerPoint());
                if (setRhoZero==0){
                    // add the error value to the density function
                    this->errorVector( gIndex) = 0.75;
                }
                else {
                    // accumulate the error value to the density function
                    this->errorVector( gIndex)  += pow(0.75, ElmAccuNmb( gIndex));//G series: sum = 3
                    ElmAccuNmb( gIndex) += 1.;
                }
            }
        }
    }
    //  End error as a function
    gsMultiPatch<> error_ml;
    error_sol.extract(error_ml);

    //------------------------------------------------------------------------
    // smoothing part of the error: we project density into more regular space
    //------------------------------------------------------------------------

    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    // Elements used for numerical integration
    A.setIntegrationElements(this->m_basis);
    A.options().setSwitch("SameElement",false);

    // Set the discretization space
    space u             = A.getSpace(this->m_basis);
    // Solution vector and solution variable
    gsMatrix<> densityVector;
    //...
    solution density_sol = A.getSolution(u, densityVector);
    //...
    //u.setup(bc_mae, dirichlet::l2Projection, 0);
    A.initSystem();

    auto rho             = A.getCoeff(error_ml);
    A.assemble(u * rho); //rhs vector
    //... density function in this case is in adaptive mesh means r o F o Psi
    densityVector        = this->Poisson.L2ProjectScalar(A.rhs());
    gsInfo << ".";

    if (!MAmapping.empty()){
    //-------------------------------------------------------------------------------------------------
    // If a given density is defined on an adaptive mesh, it must first be projected to a uniform mesh
    // since the initial mapping will be used in the composition.  (r o F o Psi) to (r o F)
    //-------------------------------------------------------------------------------------------------

    const gsKnotVector<double> kv1 =  static_cast<gsTensorNurbs<2> &>( MAmapping.patch(0)).knots(0);
    const gsKnotVector<double> kv2 =  static_cast<gsTensorNurbs<2> &>( MAmapping.patch(0)).knots(1);
    const index_t degree1 =  static_cast<gsTensorNurbs<2> &>( MAmapping.patch(0)).degree(0);
    const index_t degree2 =  static_cast<gsTensorNurbs<2> &>( MAmapping.patch(0)).degree(1);
    if (MAmapping.dim()==2){
    gsMatrix<> rhsVector = A.rhs();
    rhsVector.setZero();
    //...
    assemble_rhsvector_2d(degree1, degree2, kv1, kv2, MAmapping.patch(0).coefs(), densityVector, rhsVector);
    densityVector           = this->Poisson.L2ProjectScalar(rhsVector);
    gsInfo << ".";
    }
    else{
    const gsKnotVector<double> kv3 =  static_cast<gsTensorNurbs<3> &>( MAmapping.patch(0)).knots(2);
    const index_t degree3 =  static_cast<gsTensorNurbs<3> &>( MAmapping.patch(0)).degree(2);
    gsMatrix<> rhsVector = A.rhs();
    rhsVector.setZero();
    //...
    assemble_rhsvector_3d(degree1, degree2, degree3, kv1, kv2, kv3, MAmapping.patch(0).coefs(), densityVector, rhsVector);
    densityVector           = this->Poisson.L2ProjectScalar(rhsVector);
    gsInfo << "..";
    }
    }
    gsInfo <<densityVector.minCoeff()<<"/"<<densityVector.maxCoeff() << ".";
    //...
    gsMultiPatch<> density;
    density_sol.extract(density);
    gsInfo<<"<>\n";
    return  density;
}

/**
 * \brief Build and return a Monge–Ampère mapping as a MultiPatch object from a density function.
 *
 * \remark We solve iteratively the regularized Monge–Ampère equation
 * \f[ -\Delta \psi + 10^{-6}\psi = \left( \det(\nabla^2 \psi) - \rho(\psi) \right)^{1/d},  \f]
 * using a Picard iterative scheme.
 *
 * \param density Density function defined on the computational domain
 *                (otherwise the composition must be projected).
 * \param tolMAE Tolerance for the Picard iterative solver.
 */
void gsAdaptiveMultiPatchBuilder::buildMultiPatch(const gsMultiPatch<> &density, const double tolMAE) const
{

    // Target mapping
    gsMultiPatch<> Psi;
    gsInfo<<"<> Picard iterations";

    // initialize timing
    double slv_time(0);
    gsStopwatch timer;
    timer.restart();

    //! [Problem setup] Assembler
    gsExprAssembler<> A(1,1);
    A.options().setSwitch("SameElement",false);

    // Elements used for numerical integration
    A.setIntegrationElements(this->m_basis);
    gsExprEvaluator<> ev(A);

    // Set the geometry map
    geometryMap G    = A.getMap(identity_mp);

    // Set the discretization space
    space u          = A.getSpace(this->m_basis);

    // Set pow for BFO method dim in parameteric domain
    auto IGdim       = G.domainDim();

    // Set factor for BFO method
    auto gammaMAE    = factorial(G.domainDim());

    // Solution vector and solution variable
    gsMatrix<> solVector;
    solution u_sol   = A.getSolution(u, solVector);

    // Solution vector and solution for projecting each direction of mapping
    gsMatrix<> vsolVector;
    solution v_sol   = A.getSolution(u, vsolVector);

    // ---- manipulation of density function ----
    auto rho         = A.getCoeff(density);
    auto empldensity = (ev.max(abs(rho.val()))-ev.min(abs(rho.val())));
    double  int_uh_0 = 0.;
    double  int_uh_1 = 1.;

    if (empldensity < 1e-5|| this->m_IntensityMAE <= 1. )
    {
        gsInfo << " rho = 1.~~";
    }
    else{
        int_uh_0     = (this->m_IntensityMAE-1.)/empldensity;
        int_uh_1     = (ev.max(abs(rho.val()))-this->m_IntensityMAE*ev.min(abs(rho.val())))/empldensity;
    }
    double maxrho    = ev.max(int_uh_0*abs(rho.val()) + int_uh_1);
    // ......... End initialization for density.........

    // Initialize the system :  identity mapping as initial guess
    A.initSystem();

    // Compute the Neumann terms defined on physical space
    auto g_N         = A.getBdrFunction(G);
    
    //.. for compatibility: induced from laplace with pure Neumann boundary conditions
    auto Neumann_Int = ev.integralBdrBc(bc_mae.get("Neumann"), g_N.tr() * nv(G) );

    //... nromalisation of density function
    auto CoeffDensity= ev.integral((int_uh_0*abs(rho.val()) + int_uh_1));

    //... Compute the system matrix and right-hand side Monge-Ampere eqaution
    auto ExprMAE     = pow(pow(IGdim,IGdim)-gammaMAE+gammaMAE * CoeffDensity/(int_uh_0*abs(rho.val()) + int_uh_1), 1./IGdim);

    // compute conductivity coefficient : related to pure numann bd conditions
    auto CoeffCond   = Neumann_Int/ev.integral(ExprMAE);

    //Assemble rhs vector pof MAE-BFO method
    A.assemble(u*  CoeffCond * (-1.)*ExprMAE);

    // Compute the Neumann terms defined on physical space
    A.assembleBdr(bc_mae.get("Neumann"), u * g_N.tr() * nv(G) );

    gsInfo<< "." <<std::flush;// Assemblying done
    solVector        = A.rhs();
    solVector        = this->Poisson.solve(A.rhs());
    gsInfo<< "." <<std::flush; // Linear solving done

    if(MAmapping.empty()){

        // Initial guess for the gradient of potential function
        A.initSystem(IGdim);
        // Obtain control points for the gradient of Psi
        A.assemble(u * grad(u_sol));
        timer.restart();
        vsolVector       = Poisson.L2ProjectScalar(A.rhs().col(0));
        slv_time        += timer.stop();
        v_sol.extract(Psi);
        for(index_t Mp=0; Mp<identity_mp.dim(); ++Mp){
            gsMultiPatch<> PsiPsitp_temp;
            timer.restart();
            vsolVector       = Poisson.L2ProjectScalar(A.rhs().col(Mp));
            slv_time        += timer.stop();
            v_sol.extract(PsiPsitp_temp);
            Psi.embed(Mp+1);
            Psi.patch(0).coefs().col(Mp) = PsiPsitp_temp.patch(0).coefs().col(0);
        }
    // // ... correct boundary
    NormalProjectPts(Psi);
    }
    else
        Psi = MAmapping;

    // Picard loop
    auto  sv0 = solVector;
    solution u_lsol = A.getSolution(u, sv0);
    for(int ip{0}; ip<=m_maxIter; ++ip)
    {
        std::cout << "\r" << std::string(20,' ') << "\r<> Picard Iteration: " << ip << std::flush;
        
        //.. geometry map under update
        geometryMap PP    = A.getMap(Psi);

        //... density in new optimized mesh
        auto rho          = A.getCoeff(density, PP);

        // ... update residual
        solution u_sol    = A.getSolution(u, solVector);

        // Initialize the system
        A.initSystem();

        // Compute the system matrix and right-hand side Monge-Ampere eqaution
        auto  ExprMAE     = pow( abs(pow(div(PP).val(),IGdim) - gammaMAE*jac(PP).det()+ gammaMAE*CoeffDensity/(int_uh_0*abs(rho.val()) + int_uh_1)), 1./IGdim);

        // .. update Coeffeicient of conductivity
        auto IntegDensity = ev.integral(ExprMAE);
        CoeffCond         = Neumann_Int/IntegDensity;

        // Update rhs of MAE-BFO method
        A.assemble(u * CoeffCond * (-1.) * ExprMAE);//rhs vector

        // Compute the Neumann terms defined on physical space
        A.assembleBdr(bc_mae.get("Neumann"), u * g_N.tr() * nv(G) );

        // gsInfo<< " ." <<std::flush;// Assemblying done
        solVector         = this->Poisson.solve(A.rhs());

        // ... compute error and residual
        auto l2Res        = math::sqrt(ev.integral( ( grad(u_sol)-grad(u_lsol)).sqNorm()  ));
        sv0               = solVector;

        if ( l2Res < tolMAE || ip == m_maxIter ){
            auto L2MAE        = math::sqrt(ev.integral( pow( 1. - (int_uh_0*abs(rho.val()) + int_uh_1)*jac(PP).det()/CoeffDensity,2)  ));
            auto Ddet         = ev.min(jac(PP).det());
            // ! end Picard loop
            gsInfo  << ". L2_res  : "<<std::scientific<<l2Res
                    << ". L2_err  : "<<std::scientific<<L2MAE
                    << ". min(Jac): "<<std::fixed << std::setprecision(2)<<Ddet
                    << ". max(rho): "<<std::fixed << std::setprecision(2)<<maxrho<<".";
            break;
            } //
    A.initSystem(IGdim);
    // Obtain control points for the gradient of Psi
    A.assemble(u * grad(u_sol));
    for(index_t Mp=0; Mp<identity_mp.dim(); ++Mp){
        gsMultiPatch<> PsiPsitp_temp;
        timer.restart();
        vsolVector = Poisson.L2ProjectScalar(A.rhs().col(Mp));
        slv_time  += timer.stop();
        v_sol.extract(PsiPsitp_temp);
        Psi.patch(0).coefs().col(Mp) = PsiPsitp_temp.patch(0).coefs().col(0);
    }
    //--------------------
    //correct the boundary
    NormalProjectPts(Psi);
    // ...
    }//END for loop
    Psi.computeTopology();
    // ....
    this->MAmapping = Psi;
    timer.stop();
    gsInfo<<" CPU-time : "<<std::scientific<< slv_time   <<"<>\n";
};

/**
 * \brief Compute the L^2-projection of a composition and return it as a MultiPatch object.
 *
 * \remark The composition is projected onto a B-spline space using an L^2 projection.
 *         In 2D, a boundary correction can be applied if the approximation near the
 *         boundary is not sufficiently accurate (currently supported only in the tensor-product case).
 *
 * \param Cbasis Basis functions onto which the composition is projected.
 * \param quadValue Number of quadrature points used to evaluate the composition.
 *                  A recommended choice is degree(initial) * degree(MA mapping) + 1.
 * \param sepBoundary If true, an additional correction is applied at the boundary.
 */
gsMultiPatch<> gsAdaptiveMultiPatchBuilder::buildCompMultiPatch(const gsMultiBasis<> Cbasis, const int quadValue, const bool& sepBoundary) const 
{

    gsInfo<<"<L2> computes composition";

    // target mapping
    gsMultiPatch<> Psi;
    
    /// start timing
    double slv_time(0);
    gsStopwatch timer;
    timer.restart();

    //! [Problem setup] assembler 
    gsExprAssembler<> A(1,1);
    A.options().setReal("quA", quadValue);
    A.options().setSwitch("SameElement",false); // Very important for the composition of the two mappings
    //...
    A.setIntegrationElements(Cbasis);

    //... 
    space v        = A.getSpace(Cbasis);
    gsMatrix<> vsolVector;
    solution v_sol = A.getSolution(v, vsolVector);

    //...
    geometryMap PP = A.getMap(this->MAmapping);
    //...
    A.initSystem();
    A.assemble(v*v.tr());//Matrix for one component
    solver.compute( A.matrix() );
    // ...
    auto comp      = A.getCoeff(this->initial_mapping, PP);

    A.initSystem(this->initial_mapping.geoDim());
    A.assemble(v * comp.tr() );// blocked by this one
    vsolVector     = solver.solve(A.rhs().col(0));
    v_sol.extract(Psi);

    #pragma omp for
    for(index_t i=1; i<this->initial_mapping.geoDim(); ++i)
    {
        gsMultiPatch<> PsiVec;
        vsolVector = solver.solve(A.rhs().col(i));
        v_sol.extract(PsiVec);
        Psi.embed(i+1);
        Psi.patch(0).coefs().col(i) = PsiVec.patch(0).coefs();
    }

    // Correction step
    if (Psi.dim() == 2 && sepBoundary){
    // since the initial mapping will be used in the composition.  (F o Psi)
    const gsKnotVector<double> kv1 =  static_cast<gsTensorNurbs<2> &>( Psi.patch(0)).knots(0);
    const gsKnotVector<double> kv2 =  static_cast<gsTensorNurbs<2> &>( Psi.patch(0)).knots(1);
    const index_t degree1 =  static_cast<gsTensorNurbs<2> &>( Psi.patch(0)).degree(0);
    const index_t degree2 =  static_cast<gsTensorNurbs<2> &>( Psi.patch(0)).degree(1);
    //...
    // gsFunctionExpr<> sDir = gsFunctionExpr<>("0",2);
    // Initialize boundary conditions for Monge-Ampere solver
    // gsBoundaryConditions<> bc;
    // bc.setGeoMap(initial_mapping);
    // for ( gsMultiPatch<>::const_biterator
    //         bit = initial_mapping.bBegin(); bit != initial_mapping.bEnd(); ++bit)
    // {
    //     bc.addCondition( *bit, condition_type::dirichlet, &sDir );
    // }    
    gsSparseMatrix<> Mx  = assembleMass(Psi.basis(0).component(0));
    // eliminateDirichlet1D(boundaryConditionsForDirection(bc,0), A.options(), Mx);
    gsSparseMatrix<> My  = assembleMass(Psi.basis(0).component(1));
    // eliminateDirichlet1D(boundaryConditionsForDirection(bc,1), A.options(), My);
    // componenet by component
    #pragma omp parallel for
    for(index_t comp_nb = 0; comp_nb < Psi.geoDim(); ++comp_nb){
        // y = 0 
        gsVector<> rhsxy0(Mx.rows());
        assemble_rhsvector_1d(Psi.basis(0).component(0), degree1, kv1, 0., 1, comp_nb, rhsxy0);
        solver.compute( Mx );
        solver.setTolerance(1e-30);
        auto xsoly0 = solver.solve(rhsxy0);
        // y = 1
        gsVector<> rhsxy1(Mx.rows());
        assemble_rhsvector_1d(Psi.basis(0).component(0), degree1, kv1, 1., 1, comp_nb, rhsxy1);
        solver.compute( Mx );
        solver.setTolerance(1e-30);
        auto xsoly1 = solver.solve(rhsxy1);    
        // x = 0 
        gsVector<> rhsx0y(My.rows());
        assemble_rhsvector_1d(Psi.basis(0).component(1), degree2, kv2, 0., 0, comp_nb, rhsx0y);
        solver.compute( My );
        solver.setTolerance(1e-30);
        auto x0soly = solver.solve(rhsx0y);
        // x = 1 
        gsVector<> rhsx1y(My.rows());
        assemble_rhsvector_1d(Psi.basis(0).component(1), degree2, kv2, 1., 0, comp_nb, rhsx1y);
        solver.compute( My );
        solver.setTolerance(1e-30);
        auto x1soly = solver.solve(rhsx1y);
        CorrecBoundary(Psi, 0, comp_nb, xsoly0, xsoly1, x0soly, x1soly, true);
    }
    gsInfo << ".";
    }
    Psi.computeTopology();
    //...
    slv_time += timer.stop();
    timer.stop();
    gsInfo<<" CPU-time "<< slv_time   <<"<>\n";
    return Psi;
};

/**
 * \brief computes the projection of a composition and return a MultiPatch object :: Collocation
 *
 * \param Cbasis Basis functions onto which the composition is projected.
 */
gsMultiPatch<> gsAdaptiveMultiPatchBuilder::buildColCompMultiPatch(const gsMultiBasis<> Cbasis) const 
{

    gsInfo<<"<Col> computes composition";

    gsMultiPatch<> Psi; 
    double slv_time(0);
    gsStopwatch timer;
    timer.restart();

    gsMatrix<> initialGrid         = Cbasis.basis(0).anchors();
    // Evaluate f at the Greville points
    gsMatrix<> intervalues         = this->MAmapping.patch(0).eval(initialGrid);
    intervalues                    = intervalues.cwiseMax(0).cwiseMin(1);
    gsMatrix<> finalValues         = this->initial_mapping.patch(0).eval(intervalues);
    gsGeometry<>::uPtr interpolant = Cbasis.basis(0).interpolateData(finalValues, initialGrid);
    // extract the mapping
    Psi.addPatch(give(interpolant));
    Psi.computeTopology();
    //...
    slv_time += timer.stop();
    timer.stop();
    gsInfo<<" CPU-time "<< slv_time   <<"<>\n";
    return Psi;
};

/**
 * \brief Compute the projection (fitting) of a composition and return it as a MultiPatch object.
 *
 * \remark The composition is approximated in a B-spline space using a least-squares (L^2) fitting approach.
 *         In 2D, a boundary correction can be applied if the approximation near the boundary
 *         is not sufficiently accurate (currently supported only in the tensor-product case).
 *
 * \param Cbasis Basis functions onto which the composition is projected.
 * \param numElData Refinement level of the sampling grid used for the least-squares fitting,
 *                  relative to the identity mapping.
 * \param lambda Regularization (preconditioning) parameter; recommended value is 10^{-6},
 *               or at least smaller than 10^{-3}.
 * \param sepboundary if one want to separate boundary from inner  points in fitting procedure (based on Newton iteration( solve optimization pr))
 */
gsMultiPatch<> gsAdaptiveMultiPatchBuilder::buildFitCompMultiPatch(const gsMultiBasis<> Cbasis, const int numElData, const real_t lambda, const bool& sepboundary) const 
{

    gsInfo<<"<Fit> computes composition";
    assert(Cbasis.dim() == 2 && "Only single patch 2D fitting is implemented so far.");

    gsMultiPatch<> Psi;    
    // Copy tensor basis
    gsTHBSplineBasis<2>  THB ( Cbasis.basis(0));

    double slv_time(0);
    gsStopwatch timer;
    timer.restart();
    //...  just to generate grid in the computational domain  (we start from identity as we want tensor grid)
    gsMultiBasis<> T_tbasis(identity_mp, true);

    while ( T_tbasis.basis(0).size() <= std::max(Cbasis.basis(0).size()*30, identity_mp.basis(0).size()*numElData))
    {
        T_tbasis.uniformRefine();
    }
    gsInfo<<":gridsize="<<T_tbasis.basis(0).size()<<"/"<<Cbasis.basis(0).size();
    gsMatrix<> initialGrid    = T_tbasis.basis(0).anchors();

    // ------------------------------
    // ... generate data for fitting
    // gsMatrix<> initialGrid             = Cbasis.basis(0).anchors();
    // Evaluate f at the Greville points
    gsMatrix<> intervalues          = this->MAmapping.patch(0).eval(initialGrid);
    intervalues                     = intervalues.cwiseMax(0).cwiseMin(1);
    gsMatrix<> finalValues           = this->initial_mapping.patch(0).eval(intervalues);


    //-----------------------
    // start fitting
    //-----------------------

    // Create hierarchical refinement object
    gsFitting<> ref(initialGrid, finalValues, THB);
    //... compute coefs
    ref.compute(lambda);

    if(sepboundary){
    // number of grid points in one direction, assuming a square grid.
    index_t ngrids        = sqrt(initialGrid.cols()); 
    std::vector<index_t> boundaryIdx;
    for (index_t i = 0; i < ngrids; ++i)
    {
        // y = 0 (bottom)
        boundaryIdx.push_back(i * ngrids + 0);
        // y = 1 (top)
        boundaryIdx.push_back(i * ngrids + (ngrids - 1));
        // x = 0 (left)
        boundaryIdx.push_back(0 * ngrids + i);
        // x = 1 (right)
        boundaryIdx.push_back((ngrids - 1) * ngrids + i);
    }
    // Remove duplicates from boundary indices
    std::sort(boundaryIdx.begin(), boundaryIdx.end());
    boundaryIdx.erase(std::unique(boundaryIdx.begin(), boundaryIdx.end()),
                    boundaryIdx.end());
    // Number of interior points
    index_t totalPts = ngrids * ngrids;
    index_t nbInterior = totalPts - boundaryIdx.size();
    // Combine interior and boundary indices for fitting
    std::vector<index_t> interpIdx;
    interpIdx.push_back(nbInterior);
    for (auto b : boundaryIdx)
        interpIdx.push_back(b);
    //... initialization of geometry for fitting
    // gsMultiPatch<> PsiCol = initial_mapping;
    // while( PsiCol.basis(0).numElements() < Cbasis.basis(0).numElements())   
    //     PsiCol.uniformRefine();
    // ref.initializeGeometry(PsiCol.patch(0).coefs(), initialGrid);
    ref.parameterProjectionSepBoundary(1e-4, interpIdx);
    }

    //! [extract the mapping]
    Psi.addPatch(give(*ref.result()));
    Psi.computeTopology();
    //...
    slv_time += timer.stop();
    timer.stop();
    gsInfo<<" CPU-time "<< slv_time   <<"<>\n";
    return Psi;
};


//........................................................................
//........... useful functions for moving meshes B-spline basis ..........
//........................................................................

gsSparseMatrix<> gsAdaptiveMultiPatchBuilder::assembleMass(const gsBasis<>& basis) const
{
    gsExprAssembler<> mass(1,1);
    mass.setIntegrationDomain(basis.domain());
    typename gsExprAssembler<>::space u = mass.getSpace(basis);
    mass.initMatrix();
    mass.assemble( u * u.tr() );
    gsSparseMatrix<> result;
    mass.matrix_into(result);
    return result;
}

gsBoundaryConditions<> gsAdaptiveMultiPatchBuilder::boundaryConditionsForDirection( const gsBoundaryConditions<>& bc, index_t direction ) const
{
    gsBoundaryConditions<> result;

    for ( index_t i = 1; i <= 2; ++i)
    {
        patchSide global(0,i+2*direction), local(0,i);
        const boundary_condition<double>* cond = bc.getConditionFromSide(global);
        if (cond)
            result.addCondition(local,cond->type(),cond->function());
    }
    return result;
}

void gsAdaptiveMultiPatchBuilder::eliminateDirichlet1D(const gsBoundaryConditions<>& bc,
                          const gsOptionList& opt, gsSparseMatrix<> & result) const
{
    dirichlet::strategy ds = (dirichlet::strategy)opt.askInt("DirichletStrategy",dirichlet::elimination);
    if (ds == dirichlet::elimination)
    {
        patchSide west(0,boundary::west), east(0,boundary::east);
        index_t i = 0;
        if (bc.getConditionFromSide( west ) && bc.getConditionFromSide( west )->type() == condition_type::dirichlet ) i += 1;
        if (bc.getConditionFromSide( east ) && bc.getConditionFromSide( east )->type() == condition_type::dirichlet ) i += 2;
        if (i%2 + i/2 >= result.rows() || i%2 + i/2 >= result.cols())
            result.resize(0,0);
        else switch ( i )
             {
             case 0: break;
             case 1: result = result.block( 1, 1, result.rows()-1, result.cols()-1 ); break;
             case 2: result = result.block( 0, 0, result.rows()-1, result.cols()-1 ); break;
             case 3: result = result.block( 1, 1, result.rows()-2, result.cols()-2 ); break;
             }
    }
    else
        GISMO_ERROR("Unknown Dirichlet strategy.");
}

// Project control points following  normal direction at the boundaries for square domain (Exact square recovery after refinement)
void gsAdaptiveMultiPatchBuilder::NormalProjectPts(gsMultiPatch<>& Psi) const
{
    // normal Projection of control points (exact geometry)
    for (size_t boxNumber = 0; boxNumber < identity_mp.nPatches(); ++boxNumber)
    {
        index_t bndIter  = 1;
        #pragma omp parallel for
        for (index_t dim = 0; dim < Psi.dim(); ++dim){
            float_t lVal    = identity_mp.patch(boxNumber).coef( identity_mp.patch(boxNumber).basis().boundary(bndIter).at(0) ).array()[dim];
            // x=0 control points be like (0,:) in this case
            for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(bndIter).size(); ++i_x) 
            {
                Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(bndIter).at(i_x) ).array()[dim] = lVal;
            }
            bndIter ++;
            float_t hVal    = identity_mp.patch(boxNumber).coef( identity_mp.patch(boxNumber).basis().boundary(bndIter).at(0) ).array()[dim];
            // x=1 control points be like (1,:) in this case
            for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(bndIter).size(); ++i_x)
            {
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(bndIter).at(i_x) ).array()[dim] = hVal;
            }
            bndIter ++;
        }
    }
}

// Correct the control point at the boundary of a final mapping 
void gsAdaptiveMultiPatchBuilder::CorrecBoundary(gsMultiPatch<>& Psi, const index_t& patchNumber, const index_t& patch_cmp, const gsMatrix<>& xsoly0, const gsMatrix<>& xsoly1, const gsMatrix<>& x0soly, const gsMatrix<>& x1soly, const bool& corners) const
{
    // ...
    if (corners){
    // update corners
    for (int i_x =0; i_x < Psi.patch(patchNumber).basis().boundary(1).size(); ++i_x)
    {
    // x=0 control points be like (0,:) in this case
    Psi.patch(patchNumber).coef( Psi.patch(patchNumber).basis().boundary(1).at(i_x) ).array()[patch_cmp] = x0soly(i_x);
    // x=1 control points be like (1,:) in this case
    Psi.patch(patchNumber).coef( Psi.patch(patchNumber).basis().boundary(2).at(i_x) ).array()[patch_cmp] = x1soly(i_x);
    }

    for (int i_x =0; i_x < Psi.patch(patchNumber).basis().boundary(3).size(); ++i_x) 
    {
    // y=0 control points be like (:,0) in this case
    Psi.patch(patchNumber).coef( Psi.patch(patchNumber).basis().boundary(3).at(i_x) ).array()[patch_cmp] = xsoly0(i_x);
    // y=1 control points be like (:,1) in this case
    Psi.patch(patchNumber).coef( Psi.patch(patchNumber).basis().boundary(4).at(i_x) ).array()[patch_cmp] = xsoly1(i_x);
    }
    }
    else{
    for (int i_x =0; i_x < Psi.patch(patchNumber).basis().boundary(1).size()-1; ++i_x)
    {
    // x=0 control points be like (0,:) in this case
    Psi.patch(patchNumber).coef( Psi.patch(patchNumber).basis().boundary(1).at(i_x+1) ).array()[patch_cmp] = x0soly(i_x);
    // x=1 control points be like (1,:) in this case
    Psi.patch(patchNumber).coef( Psi.patch(patchNumber).basis().boundary(2).at(i_x+1) ).array()[patch_cmp] = x1soly(i_x);
    }

    for (int i_x =0; i_x < Psi.patch(patchNumber).basis().boundary(3).size()-1; ++i_x) 
    {
    // y=0 control points be like (:,0) in this case
    Psi.patch(patchNumber).coef( Psi.patch(patchNumber).basis().boundary(3).at(i_x+1) ).array()[patch_cmp] = xsoly0(i_x);
    // y=1 control points be like (:,1) in this case
    Psi.patch(patchNumber).coef( Psi.patch(patchNumber).basis().boundary(4).at(i_x+1) ).array()[patch_cmp] = xsoly1(i_x);
    }        
    }
}

// find the knot span index for a given x value
index_t gsAdaptiveMultiPatchBuilder::find_span(const gsKnotVector<double>& knots, const index_t& degree, const double& x)  const 
{
    index_t low = degree;
    index_t high = knots.size() - 1 - degree;
    //  Check if point is exactly on left/right boundary, or outside domain
    if (x <= knots[low]) return low;
    if (x >= knots[high]) return high - 1;
    // Binary search for the knot span index
    index_t mid = (low + high) / 2;
    while (x < knots[mid] || x >= knots[mid + 1]) {
        if (x < knots[mid]) high = mid;
        else low = mid;
        mid = (low + high) / 2;
    }
    return mid;
};

// compute the basis functions for a given knot span and degree and x value
void gsAdaptiveMultiPatchBuilder::basis_functions(const gsKnotVector<double>& knots, const index_t& degree, const double& x, index_t& span,
                     gsVector<double>& dbasis0,
                     gsVector<double>& dbasis1)  const 
{
    index_t minm = std::min(degree,2);// if degree == 1, then all second derivatives are 0

    gsVector<double> left(degree), right(degree);
    gsMatrix<double> ndu(degree + 1, degree + 1);
    gsMatrix<double> a(2, degree + 1);
    gsMatrix<double> ders(3, degree + 1);
    span     = gsAdaptiveMultiPatchBuilder::find_span(knots, degree, x);

    ndu(0,0)        = 1.0;
    // computes basis functions
    // #pragma omp parallel for
    for (index_t j = 0; j < degree; ++j) {
        left(j)      = x - knots[span - j];
        right(j)     = knots[span + j + 1] - x;
        double saved = 0.0;
        for (index_t r = 0; r <= j; ++r) {
            // compute inverse of knot differences and save them into lower triangular part of ndu
            ndu(j + 1,r)  = 1./(right(r) + left(j - r));
            //compute basis functions and save them into upper triangular part of ndu
            double temp   = ndu(r,j) * ndu(j + 1,r);
            ndu(r,j + 1)  = saved + right(r) * temp;
            saved         = left(j - r) * temp;
        }
        ndu(j + 1, j + 1) = saved;
    }
    // compute the derivatives of the basis functions
    // #pragma omp parallel for
    for (index_t r = 0; r <= degree; ++r) {
        ders(0,r) = ndu(r, degree);

        index_t s1 = 0, s2 = 1;
        a(0, 0)    = 1.0;
        for (index_t k = 1; k <= minm; ++k) {
            double d   = 0.0;
            index_t rk = r - k;
            index_t pk = degree - k;
            if (r >= k) {
                a(s2, 0) = a(s1, 0) * ndu(pk + 1, rk);
                d        = a(s2, 0) * ndu(rk, pk);
            }
            index_t j1 = (rk > -1) ? 1 : -rk;
            index_t j2 = (r - 1 <= pk) ? k - 1 : degree - r;
            for (index_t j = j1; j <= j2; ++j)
                a(s2,j) = (a(s1, j) - a(s1, j - 1)) * ndu(pk + 1, rk + j);
            for (index_t j = j1; j <= j2; ++j)
                d += a(s2, j) * ndu(rk + j, pk);
            if (r <= pk) {
                a(s2, k) = -a(s1, k - 1) * ndu(pk + 1, r);
                d += a(s2, k) * ndu(r, pk);
            }
            ders(k, r) = d;
            std::swap(s1, s2);
        }
    }
    for (index_t j = 0; j <= degree; ++j) {
        dbasis0(j) = ders(0,j);
        dbasis1(j) = ders(1,j) * degree;
    }
};

// compute the right-hand side vector for the composition of Left mapping and Monge-Ampere mapping
void gsAdaptiveMultiPatchBuilder::assemble_rhsvector_1d(
    const gsBasis<>& basis,
    const index_t& p1,
    const gsKnotVector<double>& knots_1,
    const double& valpt,
    const index_t& dir_valpt,
    const index_t& comp_nb,
    gsVector<double>& rhs) const 
{
    rhs.setZero();

    const double pi = 3.141592653589793;
    index_t m = 8*p1 + 1;
    index_t nRoots = (m + 1) / 2;

    static gsVector<double> u, w;
    static bool initialized = false;

    if (!initialized)
    {
        u.resize(m);
        w.resize(m);

        for (index_t i = 0; i < nRoots; ++i)
        {
            double t = std::cos(pi*(i + 0.75)/(m + 0.5));

            for (index_t j = 0; j < 30; ++j)
            {
                double p0 = 1.0, p1_ = t;

                for (index_t k = 1; k < m; ++k)
                {
                    double pn = ((2.0*k+1.0)*t*p1_ - k*p0)/(k+1.0);
                    p0 = p1_;
                    p1_ = pn;
                }

                double denom = (1.0 - t*t);
                if (std::abs(denom) < 1e-14) break;

                double dp = m*(p0 - t*p1_)/denom;
                double dt = -p1_/dp;

                t += dt;

                if (std::abs(dt) < 1e-14)
                {
                    u(i) = t;
                    u(m-i-1) = -t;
                    w(i) = 2.0/((1.0 - t*t)*dp*dp);
                    w(m-i-1) = w(i);
                    break;
                }
            }
        }
        initialized = true;
    }

    index_t nb1 = knots_1.size() - p1 - 1;
    index_t ne1 = nb1 - p1;

    gsVector<double> points_cmp(2);
    points_cmp(dir_valpt) = valpt;
    index_t dir_op = 1 - dir_valpt;

    for (index_t ie1 = 0; ie1 < ne1; ++ie1)
    {
        double a1 = knots_1[ie1 + p1];
        double b1 = knots_1[ie1 + p1 + 1];

        if (std::abs(b1 - a1) < 1e-14)
            continue;

        double c0 = 0.5*(a1 + b1);
        double c1 = 0.5*(b1 - a1);

        for (index_t g1 = 0; g1 < m; ++g1)
        {
            double x1 = c1*u(g1) + c0;
            double w1 = c1*w(g1);

            index_t span1;
            gsVector<double> N(p1+1), dN(p1+1);

            basis_functions(knots_1, p1, x1, span1, N, dN);

            points_cmp(dir_op) = x1;
            // gsInfo << "before :" << points_cmp;
            gsMatrix<> val_i = MAmapping.patch(0).eval(points_cmp);
            val_i(dir_valpt) = valpt;
            val_i.cwiseMax(0).cwiseMin(1);
            // gsInfo << "after :" << val_i << "---\n";
            double val_f = initial_mapping.patch(0).eval(val_i)[comp_nb];

            for (index_t i = 0; i <= p1; ++i)
            {
                index_t gi = span1 - p1 + i;
                rhs[gi] += val_f * N(i) * w1;
            }
        }
    }
};
// End of the function

// compute the right-hand side vector for the adaptive multi-patch assembly : one patch
// This function computes the right-hand side vector for a given knot span and degree
// vector_un : vector of control points for the solution
// vector_mp : vector of control points for the Monge-Ampere mapping : unit square
void gsAdaptiveMultiPatchBuilder::assemble_rhsvector_2d(const index_t& p1, const index_t& p2,
                           const gsKnotVector<double>& knots_1, const gsKnotVector<double>& knots_2,
                           const gsMatrix<double>& vector_mp, const gsMatrix<double>& vector_un,
                           gsMatrix<double>& rhs) const {
    
    rhs.setZero();
    const double pi = 3.141592653589793;
    index_t m       = p1 + p2 + 1;
    index_t nRoots  = (m + 1) / 2;

    gsVector<double> w(m);
    gsVector<double> u(m);

    #pragma omp parallel for
    for (index_t i = 0; i < nRoots; ++i) {
        double t = std::cos(pi* (i + 0.75) / (m + 0.5));
        for (index_t j = 0; j < 30; ++j) {
            double p_0   = 1.0;
            double p_1   = t;
            for (index_t k = 1; k < m; ++k) {
                double pn = ((2.0 * k + 1.0) * t * p_1 - k * p_0) / (1.0 + k);
                p_0       = p_1;
                p_1       = pn;
            }
            double dp = m * (p_0 - t * p_1) / (1.0 - t * t);
            double dt = -p_1 / dp;
            t += dt;
            if (std::abs(dt) < 1e-14) {
                u(i) = t;
                u(m - i - 1) = -t;
                w(i) = 2.0 / ((1.0 - t * t) * dp * dp);
                w(m - i - 1) = w(i);
                break;
            }
        }
    }
    // Compute the number of basis functions in each direction
    index_t nb1   = knots_1.size() - p1 - 1;
    index_t nb2   = knots_2.size() - p2 - 1;

    // Compute the number of elements in each direction
    index_t ne1 = nb1 - p1;
    index_t ne2 = nb2 - p2;

    #pragma omp parallel for
    for (index_t ie1 = 0; ie1 < ne1; ++ie1) {

        // coefs of the element for quadrature points
        double a1   = knots_1[ie1 + p1];
        double b1   = knots_1[ie1 + p1 + 1];
        double c0_1 = 0.5 * (a1 + b1);
        double c1_1 = 0.5 * (b1 - a1);

        for (index_t ie2 = 0; ie2 < ne2; ++ie2) {

            double a2   = knots_2[ie2 + p2];
            double b2   = knots_2[ie2 + p2 + 1];
            double c0_2 = 0.5 * (a2 + b2);
            double c1_2 = 0.5 * (b2 - a2);

            for (index_t g1 = 0; g1 < m; ++g1) {
                // map the quadrature points to the element
                double x1 = c1_1 * u(g1) + c0_1;
                double w1 = c1_1 * w(g1);
                // Compute the basis functions in the first direction
                index_t span1;
                gsVector<double> xbasis_0(p1 + 1);
                gsVector<double> xbasis_1(p1 + 1);                
                basis_functions(knots_1, p1, x1, span1, xbasis_0, xbasis_1);

                for (index_t g2 = 0; g2 < m; ++g2) {
                    // map the quadrature points to the element
                    double x2 = c1_2 * u(g2) + c0_2;
                    double w2 = c1_2 * w(g2);
                    // Compute the basis functions in the second direction
                    index_t span2;
                    gsVector<double> ybasis_0(p2 + 1);
                    gsVector<double> ybasis_1(p2 + 1);
                    basis_functions(knots_2, p2, x2, span2, ybasis_0, ybasis_1);

                    // Assembles solution in uniform mesh
                    double val_un = 0.0;
                    double ad_x   = 0.0;
                    double ad_xx  = 0.0;
                    double ad_xy  = 0.0;
                    double ad_y   = 0.0;
                    double ad_yx  = 0.0;
                    double ad_yy  = 0.0;
                    for (index_t i = 0; i <= p1; ++i) {
                        for (index_t j = 0; j <= p2; ++j) {
                            index_t gi  = (span1 - p1+i) + nb1*(span2 - p2+j);

                            double bi_0 = xbasis_0(i) * ybasis_0(j);

                            double bi_x = xbasis_1(i) * ybasis_0(j);
                            double bj_y = xbasis_0(i) * ybasis_1(j);

                            ad_x       += vector_mp.col(0)(gi) * bi_0;
                            ad_y       += vector_mp.col(1)(gi) * bi_0;

                            ad_xx      += vector_mp.col(0)(gi) * bi_x;
                            ad_xy      += vector_mp.col(0)(gi) * bj_y;
                            ad_yx      += vector_mp.col(1)(gi) * bi_x;
                            ad_yy      += vector_mp.col(1)(gi) * bj_y;

                            val_un     += vector_un(gi) * bi_0;
                        }
                    }
                    // basis functions in the image of quadrature points by Monge-Ampere mapping
                    index_t ad_span1;
                    gsVector<double> ad_xbasis_0(p1 + 1);
                    gsVector<double> ad_xbasis_1(p1 + 1);
                    basis_functions(knots_1, p1, ad_x, ad_span1, ad_xbasis_0, ad_xbasis_1);
                    index_t ad_span2;
                    gsVector<double> ad_ybasis_0(p2 + 1);
                    gsVector<double> ad_ybasis_1(p2 + 1);
                    basis_functions(knots_2, p2, ad_y, ad_span2, ad_ybasis_0, ad_ybasis_1);

                    double weight = w1 * w2 * std::abs(ad_xx * ad_yy - ad_xy * ad_yx);
                    for (index_t i = 0; i <= p1; ++i) {
                        for (index_t j = 0; j <= p2; ++j) {
                            index_t gi   = (ad_span1 - p1+i )+ nb1*(ad_span2 - p2+j);
                            double bi_0  = ad_xbasis_0(i) * ad_ybasis_0(j);

                            rhs(gi) += val_un * bi_0 * weight;
                        }
                    }
                }
            }
        }
    }
};
// End of the function

// compute the right-hand side vector for the adaptive multi-patch assembly : one patch
// This function computes the right-hand side vector for a given knot span and degree
// vector_un : vector of control points for the solution
// vector_mp : vector of control points for the Monge-Ampere mapping : unit square
void gsAdaptiveMultiPatchBuilder::assemble_rhsvector_3d(const index_t& p1, const index_t& p2, const index_t& p3,
                           const gsKnotVector<double>& knots_1, const gsKnotVector<double>& knots_2, const gsKnotVector<double>& knots_3,
                           const gsMatrix<double>& vector_mp, const gsMatrix<double>& vector_un,
                           gsMatrix<double>& rhs) const {
    
    rhs.setZero();
    const double pi = 3.141592653589793;
    index_t m       = p1 + p2 + 1;
    index_t nRoots  = (m + 1) / 2;

    gsVector<double> w(m);
    gsVector<double> u(m);

    #pragma omp parallel for
    for (index_t i = 0; i < nRoots; ++i) {
        double t = std::cos(pi* (i + 0.75) / (m + 0.5));
        for (index_t j = 0; j < 30; ++j) {
            double p_0   = 1.0;
            double p_1   = t;
            for (index_t k = 1; k < m; ++k) {
                double pn = ((2.0 * k + 1.0) * t * p_1 - k * p_0) / (1.0 + k);
                p_0       = p_1;
                p_1       = pn;
            }
            double dp = m * (p_0 - t * p_1) / (1.0 - t * t);
            double dt = -p_1 / dp;
            t += dt;
            if (std::abs(dt) < 1e-14) {
                u(i) = t;
                u(m - i - 1) = -t;
                w(i) = 2.0 / ((1.0 - t * t) * dp * dp);
                w(m - i - 1) = w(i);
                break;
            }
        }
    }
    // Compute the number of basis functions in each direction
    index_t nb1 = knots_1.size() - p1 - 1;
    index_t nb2 = knots_2.size() - p2 - 1;
    index_t nb3 = knots_3.size() - p3 - 1;
    index_t nb12 = nb1 * nb2;
    // Compute the number of elements in each direction
    index_t ne1 = nb1 - p1;
    index_t ne2 = nb2 - p2;
    index_t ne3 = nb3 - p3;

    #pragma omp parallel for
    for (index_t ie1 = 0; ie1 < ne1; ++ie1) {

        // coefs of the element for quadrature points
        double a1   = knots_1[ie1 + p1];
        double b1   = knots_1[ie1 + p1 + 1];
        double c0_1 = 0.5 * (a1 + b1);
        double c1_1 = 0.5 * (b1 - a1);

        for (index_t ie2 = 0; ie2 < ne2; ++ie2) {

            double a2   = knots_2[ie2 + p2];
            double b2   = knots_2[ie2 + p2 + 1];
            double c0_2 = 0.5 * (a2 + b2);
            double c1_2 = 0.5 * (b2 - a2);

            for (index_t ie3 = 0; ie3 < ne3; ++ie3) {

                double a3   = knots_3[ie3 + p3];
                double b3   = knots_3[ie3 + p3 + 1];
                double c0_3 = 0.5 * (a3 + b3);
                double c1_3 = 0.5 * (b3 - a3);

                for (index_t g1 = 0; g1 < m; ++g1) {
                    // map the quadrature points to the element
                    double x1 = c1_1 * u(g1) + c0_1;
                    double w1 = c1_1 * w(g1);
                    // Compute the basis functions in the first direction
                    index_t span1;
                    gsVector<double> xbasis_0(p1 + 1);
                    gsVector<double> xbasis_1(p1 + 1);                
                    basis_functions(knots_1, p1, x1, span1, xbasis_0, xbasis_1);

                    for (index_t g2 = 0; g2 < m; ++g2) {
                        // map the quadrature points to the element
                        double x2 = c1_2 * u(g2) + c0_2;
                        double w2 = c1_2 * w(g2);
                        // Compute the basis functions in the second direction
                        index_t span2;
                        gsVector<double> ybasis_0(p2 + 1);
                        gsVector<double> ybasis_1(p2 + 1);
                        basis_functions(knots_2, p2, x2, span2, ybasis_0, ybasis_1);
                        for (index_t g3 = 0; g3 < m; ++g3) {
                            // map the quadrature points to the element
                            double x3 = c1_3 * u(g3) + c0_3;
                            double w3 = c1_3 * w(g3);
                            // Compute the basis functions in the second direction
                            index_t span3;
                            gsVector<double> zbasis_0(p3 + 1);
                            gsVector<double> zbasis_1(p3 + 1);
                            basis_functions(knots_3, p3, x3, span3, zbasis_0, zbasis_1);

                            // Assembles solution in uniform mesh // reshape
                            double val_un = 0.0;
                            double ad_x   = 0.0;
                            double ad_xx  = 0.0;
                            double ad_xy  = 0.0;
                            double ad_xz  = 0.0;
                            double ad_y   = 0.0;
                            double ad_yx  = 0.0;
                            double ad_yy  = 0.0;
                            double ad_yz  = 0.0;
                            double ad_z   = 0.0;
                            double ad_zx  = 0.0;
                            double ad_zy  = 0.0;
                            double ad_zz  = 0.0;
                            for (index_t i = 0; i <= p1; ++i) {
                                for (index_t j = 0; j <= p2; ++j) {
                                    for (index_t k = 0; k <= p3; ++k) {
                                        index_t gi  = (span1 - p1+i) + nb1*(span2 - p2+j) + nb12*(span3 - p3+k);

                                        double bi_0 = xbasis_0(i) * ybasis_0(j) * zbasis_0(k);

                                        double bi_x = xbasis_1(i) * ybasis_0(j) * zbasis_0(k);
                                        double bj_y = xbasis_0(i) * ybasis_1(j) * zbasis_0(k);
                                        double bj_z = xbasis_0(i) * ybasis_0(j) * zbasis_1(k);

                                        ad_x       += vector_mp.col(0)(gi) * bi_0;
                                        ad_y       += vector_mp.col(1)(gi) * bi_0;
                                        ad_z       += vector_mp.col(2)(gi) * bi_0;

                                        ad_xx      += vector_mp.col(0)(gi) * bi_x;
                                        ad_xy      += vector_mp.col(0)(gi) * bj_y;
                                        ad_xz      += vector_mp.col(0)(gi) * bj_z;

                                        ad_yx      += vector_mp.col(1)(gi) * bi_x;
                                        ad_yy      += vector_mp.col(1)(gi) * bj_y;
                                        ad_yz      += vector_mp.col(1)(gi) * bj_z;

                                        ad_zx      += vector_mp.col(2)(gi) * bi_x;
                                        ad_zy      += vector_mp.col(2)(gi) * bj_y;
                                        ad_zz      += vector_mp.col(2)(gi) * bj_z;

                                        val_un     += vector_un(gi) * bi_0;
                                    }
                                }
                            }
                            // basis functions in the image of quadrature points by Monge-Ampere mapping
                            index_t ad_span1;
                            gsVector<double> ad_xbasis_0(p1 + 1);
                            gsVector<double> ad_xbasis_1(p1 + 1);
                            basis_functions(knots_1, p1, ad_x, ad_span1, ad_xbasis_0, ad_xbasis_1);
                            index_t ad_span2;
                            gsVector<double> ad_ybasis_0(p2 + 1);
                            gsVector<double> ad_ybasis_1(p2 + 1);
                            basis_functions(knots_2, p2, ad_y, ad_span2, ad_ybasis_0, ad_ybasis_1);
                            index_t ad_span3;
                            gsVector<double> ad_zbasis_0(p3 + 1);
                            gsVector<double> ad_zbasis_1(p3 + 1);
                            basis_functions(knots_3, p3, ad_z, ad_span3, ad_zbasis_0, ad_zbasis_1);

                            double weight = w1 * w2 * w3 * std::abs(ad_xx*(ad_yy*ad_zz-ad_zy*ad_yz) - ad_xy*(ad_xx*ad_zz - ad_zx*ad_yz) +ad_xz*(ad_yx*ad_zy -ad_zx*ad_yy));
                            for (index_t i = 0; i <= p1; ++i) {
                                for (index_t j = 0; j <= p2; ++j) {
                                    for (index_t k = 0; k <= p3; ++k) {
                                        index_t gi   = (ad_span1 - p1+i )+ nb1*(ad_span2 - p2+j) + nb12*(ad_span3 - p3+k);
                                        double bi_0  = ad_xbasis_0(i) * ad_ybasis_0(j) * ad_zbasis_0(k);

                                        rhs(gi) += val_un * bi_0 * weight;
                                    }
                                }
                            }
                        }//for g3
                    }//for g2
                }//for g1
            }//for ne3
        }//for ne2
    }//for ne1
};
// End of the function

// /**
//  * \brief  computes the projection of a the inverse mapping and return a MultiPatch object :: fitting
//  * 
//  * \remark The solution is approximated in a B-spline space using a least-squares (L^2) fitting approach.
//  *         solution(tn) = soFo\Psi^[n} -->  solution(tn) = soFo\Psi^{n+1}
//  *
//  * \param lastMAEmapping mapping at time step tn
//  * \param numElData Refinement level of the sampling grid used for the least-squares fitting,
//  *                  relative to the identity mapping.
//  * \param lambda Regularization (preconditioning) parameter; recommended value is 10^{-6},
//  *               or at least smaller than 10^{-3}.
//  * \param UpdateInTime Computes the composition of Psi_n^{-1} o Psi_{n+1}
//  * \todo not working yet 
//  */
// gsMultiPatch<> gsAdaptiveMultiPatchBuilder::buildInverseMultiPatch(const gsMultiPatch<> lastMAEmapping, const int numElData, const real_t lambda, const bool UpdateInTime) const 
// {

//     gsInfo<<"<Fit> computes composition";
//     gsMultiPatch<> PsiInv;
//     // Copy tensor basis
//     gsTHBSplineBasis<2>  THB ( MAmapping.basis(0));

//     double slv_time(0);
//     gsStopwatch timer;
//     timer.restart();
//     //...  just to generate grid
//     gsMultiBasis<> T_tbasis(identity_mp, true);

//     while ( T_tbasis.basis(0).numElements() <= MAmapping.basis(0).numElements()*numElData)
//     {
//         T_tbasis.uniformRefine();
//     }
//     gsInfo<<":gridsize="<<T_tbasis.basis(0).numElements()<<"/"<<MAmapping.basis(0).numElements();
//     // Inverse mappinbg first 
//     gsMatrix<> initialGrid     = T_tbasis.basis(0).anchors();
//     // Evaluate f at the Greville points
//     gsMatrix<> finalValues     = lastMAEmapping.patch(0).eval(initialGrid);
//     finalValues                = finalValues.cwiseMax(0).cwiseMin(1);

//     //! [Create  Hfitter]
//     // Create hierarchical refinement object
//     gsFitting<> ref(finalValues,  initialGrid, THB);
//     //... compute coefs
//     ref.compute(lambda);
     
//     //! [extract the mapping]
//     PsiInv.addPatch(give(*ref.result()));
//     PsiInv.computeTopology();
//     //...
//     slv_time += timer.stop();
//     timer.stop();
//     if (UpdateInTime){
//         //Computes the composition of Psi_n^{-1} o Psi_{n+1}
//         gsMultiPatch<> PsiComp;
//         // Evaluate f at the Greville points
//         gsMatrix<> IValues     = this->MAmapping.patch(0).eval(initialGrid);// Psi_{n+1}
//         IValues                = IValues.cwiseMax(0).cwiseMin(1);
//         gsMatrix<> finalValues = PsiInv.patch(0).eval(IValues);// Psi_n^{-1}
//         finalValues            = finalValues.cwiseMax(0).cwiseMin(1);
//         //! [Create  Hfitter]
//         // Create hierarchical refinement object
//         gsFitting<> ref(initialGrid, finalValues, THB);
//         //... compute coefs
//         ref.compute(lambda);
        
//         //! [extract the mapping]
//         PsiComp.addPatch(give(*ref.result()));
//         PsiComp.computeTopology();        
//         gsInfo<<" CPU-time "<< slv_time   <<"<cmp>\n";
//         return PsiComp;
//     }
//     else{
//         gsInfo<<" CPU-time "<< slv_time   <<"<>\n";
//         return PsiInv;
//     }
// };