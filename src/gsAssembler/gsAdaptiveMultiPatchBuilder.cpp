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

namespace {
// Transfer the density to the uniform mesh using
// rhs_i = integral rho(xi) B_i(Psi(xi)) |det(DPsi(xi))| dxi.
template <short_t d>
void assembleDensityRhs(const gsBasis<>& basis, const gsMatrix<>& mappingCoefs,
                        const gsMatrix<>& densityCoefs, gsMatrix<>& rhs)
{
    GISMO_ENSURE(basis.dim() == d && mappingCoefs.rows() == basis.size()
                 && mappingCoefs.cols() == d && densityCoefs.rows() == basis.size()
                 && densityCoefs.cols() == 1,
                 "Mapping and scalar density must use the same basis.");
    rhs.setZero(basis.size(), 1);

    // Use the library quadrature and basis evaluators in both 2D and 3D.
    gsGaussRule<> rule(basis, 2., 1);
    gsFuncData<> source(NEED_VALUE | NEED_DERIV | NEED_ACTIVE | SAME_ELEMENT);
    gsFuncData<> target(NEED_VALUE | NEED_ACTIVE);
    gsMatrix<> points, coefficients, values;
    gsVector<> weights;
    gsMatrix<real_t, d, d> jacobian;

    const auto domain = basis.domain();
    for (auto it = domain->beginAll(); it < domain->endAll(); ++it)
    {
        rule.mapTo(it.lowerCorner(), it.upperCorner(), points, weights);

        // Cache source values and derivatives once for the whole element.
        // The mapping and density share these basis values and active indices.
        basis.compute(points, source);
        const index_t numActive = source.actives.rows();
        coefficients.resize(numActive, d + 1);
        for (index_t a = 0; a < numActive; ++a)
        {
            const index_t i = source.actives(a, 0);
            coefficients.row(a).head(d) = mappingCoefs.row(i);
            coefficients(a, d) = densityCoefs(i);
        }
        values.noalias() = coefficients.transpose() * source.values[0];

        // Mapped points can cross knot spans, so do not use SAME_ELEMENT here.
        // Only values are needed at the image points, not derivatives.
        basis.compute(values.topRows(d), target);
        for (index_t q = 0; q < points.cols(); ++q)
        {
            jacobian.noalias() = coefficients.leftCols(d).transpose()
                               * source.deriv(q).transpose();
            const real_t weight = weights(q) * math::abs(jacobian.determinant())
                                * values(d, q);

            // Scatter serially: different source elements may update the same dof.
            for (index_t a = 0; a < target.actives.rows(); ++a)
                rhs(target.actives(a, q)) += weight * target.values[0](a, q);
        }
    }
}

// Density restricted to one side of the unit square, as a function of the tangential parameter
class EdgeFunction : public gsFunction<real_t>
{
public:
    EdgeFunction(const gsFunction<real_t>& rho, index_t dir, real_t fixedVal)
    : m_rho(rho), m_dir(dir), m_fixed(fixedVal) { }

    short_t domainDim() const { return 1; }
    short_t targetDim() const { return 1; }

    void eval_into(const gsMatrix<real_t>& u, gsMatrix<real_t>& result) const
    {
        gsMatrix<real_t> pts(2, u.cols());
        pts.row(m_dir)   = u.row(0);
        pts.row(1-m_dir).setConstant(m_fixed);
        m_rho.eval_into(pts, result);
    }
private:
    const gsFunction<real_t>& m_rho;
    index_t m_dir;
    real_t  m_fixed;
};
}


//+++++++++++++++++++++++++++
// Constructor implementation
gsAdaptiveMultiPatchBuilder::gsAdaptiveMultiPatchBuilder(const gsMultiPatch<> mapping,
                            index_t numRefine,
                            index_t maxIter,
                            double IntensityMAE,
                            index_t numReduce,
                            index_t numElevate,
                            bool exactGeo)
{
    gsInfo<<"\n <>r-refinement ";
    gsInfo<<"\n";
    gsInfo<<"\n";
    gsInfo << "/*/*/*/*/*/*/*/*/*/*/Important:::: Use --exactGeo to keep the exact boundary!\n";
    gsInfo<<"for more accuracy refine: DoFs > 10^3\n";
    gsInfo<<"\n";
    // Store the initial mapping
    this->initial_mapping= mapping; 
    // basis function for the composition
    gsMultiBasis<> mapping_basis(mapping, false);
    //-------------------------------------------------------------------------------------
    // MAE needs tensor B-splines; for a hierarchy, start from its level-0 grid.
    gsMultiBasis<> dbasis;
    for (size_t p = 0; p < mapping.nPatches(); ++p)
    {
        const gsBasis<> & basis = mapping.basis(p).source(); // Remove rational weights.
        if (const auto * h = dynamic_cast<const gsHTensorBasis<2> *>(&basis))
            dbasis.addBasis(h->tensorLevel(0).clone());
        else if (const auto * h = dynamic_cast<const gsHTensorBasis<3> *>(&basis))
            dbasis.addBasis(h->tensorLevel(0).clone());
        else
        {
            GISMO_ENSURE(dynamic_cast<const gsTensorBSplineBasis<2> *>(&basis)
                         || dynamic_cast<const gsTensorBSplineBasis<3> *>(&basis),
                         "MAE requires a tensor B-spline or hierarchical tensor basis in 2D or 3D.");
            dbasis.addBasis(basis.clone());
        }
    }
    dbasis.setTopology(mapping.topology());
    //... refine basis for convergence
    for(index_t i = 0; i < numRefine; ++i)
    {   
        dbasis.uniformRefine();
        mapping_basis.uniformRefine();
    }
    // Reduce degree if possible while maintaining minimum degree of 1
    if (numReduce > 0){
        int reduceDegree = std::min(dbasis.degree()-1, numReduce);
        dbasis.degreeDecrease(reduceDegree);
    }
    else if (numElevate > 0){
        // Elevate degree if possible
        dbasis.degreeElevate(numElevate);
    }
    // All patches share one topology-free basis on the unit square (the geometry may glue sides, e.g. closed surfaces)
    if (mapping.nPatches() > 1 || exactGeo)
        dbasis = gsMultiBasis<>(dbasis.basis(0));

    // Store input parameters
    this->m_basis        = dbasis;
    this->mapping_basis  = mapping_basis;
    this->m_maxIter      = maxIter;
    this->m_IntensityMAE = IntensityMAE;
    this->m_exactGeo     = exactGeo;
    this->DoFs           = m_basis.size() * mapping.nPatches();

    gsInfo << "nb patches " << mapping.nPatches() << " Using B-splines of degree " << dbasis.degree() << " DoFs ";

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
    this->PoissonDir     = dirichletPoissonSolver();
    gsInfo<<this->DoFs <<"<> \n";
}


// Fast diagonalization solver for the Poisson problem with homogeneous Dirichlet conditions (multipatch case);
// the Neumann version used for Monge-Ampere is `Poisson`.
gsPatchPreconditionersCreator<double>::Poisson_FastDiag gsAdaptiveMultiPatchBuilder::dirichletPoissonSolver() const
{
    gsFunctionExpr<> zero("0.", identity_mp.domainDim());
    gsBoundaryConditions<> bc;
    bc.setGeoMap(identity_mp);
    for (gsMultiPatch<>::const_biterator bit = identity_mp.bBegin(); bit != identity_mp.bEnd(); ++bit)
        bc.addCondition(*bit, condition_type::dirichlet, &zero);
    gsExprAssembler<> A(1,1);
    return gsPatchPreconditionersCreator<double>::Poisson_FastDiag(this->m_basis.basis(0), bc, A.options(), 0.);
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
    this->DoFs      = m_basis.size() * initial_mapping.nPatches();
    
    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    
    //::::::::::::::::::::      Poisson fast diagonalization solver         :::::::::::::::::::::::::
    gsPatchPreconditionersCreator<double>::Poisson_FastDiag Poisson(this->m_basis.basis(0), bc_mae, A.options(), 1e-6);  
    this->Poisson   = Poisson;
    this->PoissonDir = dirichletPoissonSolver();
    gsInfo << "<>"<< this->DoFs <<" DoFs after uniRefine <>\n";
}


/** \brief Build and return a density as a MultiPatch object from analytical function 
 * \remark (we avoid three compositions (r o F o Psi) here to be r o Psi)
 * \remark In the multipatch case one density patch (on the unit square) is built per geometry patch.
 * \param f density function defined on the physical domaine.
 */
gsMultiPatch<> gsAdaptiveMultiPatchBuilder::buildAnalyticDensity(const gsFunctionExpr<>   &f) const 
{
    gsInfo<<"<>density function: ";

    if (this->initial_mapping.nPatches() > 1)
    {
        // Global L2 projection with C0 coupling at the interfaces: the density is continuous across patches
        gsMultiBasis<>::BasisContainer bases;
        for (size_t n = 0; n < this->initial_mapping.nPatches(); ++n)
            bases.push_back(this->m_basis.basis(0).clone().release());
        gsMultiBasis<> mb(bases, this->initial_mapping);

        gsExprAssembler<> Ag(1,1);
        Ag.options().setSwitch("SameElement",false);
        Ag.setIntegrationElements(mb);
        space ug = Ag.getSpace(mb);
        gsBoundaryConditions<> bcEmpty;
        ug.setup(bcEmpty, dirichlet::l2Projection, 0);
        geometryMap Gg = Ag.getMap(this->initial_mapping);
        auto fg        = Ag.getCoeff(f, Gg);
        Ag.initSystem();
        Ag.assemble(ug * ug.tr(), ug * fg.val());

        gsMatrix<> densityVec = gsSparseSolver<>::CGDiagonal().compute(Ag.matrix()).solve(Ag.rhs());
        solution density_sol  = Ag.getSolution(ug, densityVec);
        gsMultiPatch<> density;
        density_sol.extract(density);
        gsInfo<< densityVec.maxCoeff()<< "/"<<densityVec.minCoeff()<<" <>\n";
        return density;
    }

    //! Assemnbler
    gsExprAssembler<> A(1,1);

    // Elements used for numerical integration
    A.setIntegrationElements(this->m_basis);

    //.. Evaluator
    gsExprEvaluator<> ev(A);

    // Set the discretization space
    space u             = A.getSpace(this->m_basis);

    // Solution vector and solution variable
    gsMatrix<> densityVector;
    solution density_sol = A.getSolution(u, densityVector);
    gsMultiPatch<> density;

    for (size_t n = 0; n < this->initial_mapping.nPatches(); ++n)
    {
        // Set the Target geometry map
        geometryMap GLeft   = A.getMap(this->initial_mapping.patch(n));

        // Set the source term with respect to target geometry
        auto ff             = A.getCoeff(f, GLeft);

        //u.setup(bc_mae, dirichlet::l2Projection, 0);
        A.initSystem();
        A.assemble(u* ff.val()); //rhs vector

        // prjection L2 of the composition
        densityVector        = this->Poisson.L2ProjectScalar(A.rhs());

        // save density funcrtion as multipatch
        gsMultiPatch<> density_n;
        density_sol.extract(density_n);
        density.addPatch(density_n.patch(0));
        gsInfo<< densityVector.maxCoeff()<< "/"<<densityVector.minCoeff()<<" ";
    }
    gsInfo<<"<>\n";
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
    const size_t nP = this->initial_mapping.nPatches();
    GISMO_ENSURE(nP == 1 || identity_mp.domainDim() == 2, "Multipatch density is implemented in 2D only.");
    // one unit-square basis per patch, not coupled across interfaces (piecewise constants are discontinuous)
    gsMultiBasis<>::BasisContainer bases0;
    for (size_t pn=0; pn < nP; ++pn)
        bases0.push_back(identity_mp.basis(0).clone().release());
    gsMultiBasis<> basis_0(bases0, gsBoxTopology(identity_mp.domainDim(), nP));
    basis_0.uniformRefine(setRhogrid);
    while (basis_0.basis(0).size() <= std::min<index_t>( this->m_basis.basis(0).size(), Hbasis.size()/static_cast<index_t>(nP))) // refine until having enough resolution for error representation
        basis_0.uniformRefine();

    // ... We want each element to be reprensted by one basis for all patches
    for (size_t pn=0; pn < nP; ++pn ) 
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
    const index_t nDofs0 = basis_0.basis(0).size();
    for (size_t pn=0; pn < nP; ++pn )// for all patches
    {
        // for all elements in patch pn
        typename gsBasis<>::domainIter domIt =  // add patchInd to domainiter ?
            Hbasis.basis(pn).domain()->beginAll();
        typename gsBasis<>::domainIter domItEnd =  // add patchInd to domainiter ?
            Hbasis.basis(pn).domain()->endAll();
        for (; domIt<domItEnd; ++domIt )
        {
            if( elMarked[ globalCount++ ] ){ // refine this element ?
                // dof index in the basis_0 (patches are stored one after another)
                auto gIndex = static_cast<index_t>(pn)*nDofs0 + basis_0.basis(pn).elementIndex(domIt.centerPoint());
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
    solution density_sol = A.getSolution(u, densityVector);
    gsMultiPatch<> density;

    for (size_t pn=0; pn < nP; ++pn)
    {
    A.initSystem();

    auto rho             = A.getCoeff(error_ml.patch(pn));
    A.assemble(u * rho); //rhs vector
    //... density function in this case is in adaptive mesh means r o F o Psi
    densityVector        = this->Poisson.L2ProjectScalar(A.rhs());
    gsInfo << ".";

    if (!MAmapping.empty()){
    //-------------------------------------------------------------------------------------------------
    // If a given density is defined on an adaptive mesh, it must first be projected to a uniform mesh
    // since the initial mapping will be used in the composition.  (r o F o Psi) to (r o F)
    //-------------------------------------------------------------------------------------------------

    // Use the actual shared basis; the MAE/MMPDE mapping need not be a NURBS.
    const gsBasis<>& basis = MAmapping.patch(pn).basis();
    gsMatrix<> rhsVector;
    if (basis.dim() == 2)
        assemble_rhsvector_2d(basis, MAmapping.patch(pn).coefs(), densityVector, rhsVector);
    else
        assemble_rhsvector_3d(basis, MAmapping.patch(pn).coefs(), densityVector, rhsVector);
    densityVector = this->Poisson.L2ProjectScalar(rhsVector);
    gsInfo << ".";
    }
    gsInfo <<densityVector.minCoeff()<<"/"<<densityVector.maxCoeff() << ".";
    gsMultiPatch<> density_n;
    density_sol.extract(density_n);
    density.addPatch(density_n.patch(0));
    }
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
 */
void gsAdaptiveMultiPatchBuilder::buildMultiPatch(const gsMultiPatch<> &density) const
{
    GISMO_ENSURE(m_maxIter >= 0, "Picard update budget must be nonnegative.");
    if (this->m_exactGeo || this->initial_mapping.nPatches() > 1)
    {
        buildMultiPatchMMPDE(density);
        return;
    }

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
    setBoundaryControlPointsAlongNormal(Psi);
    }
    else
        Psi = MAmapping;

    // Picard loop
    auto  sv0 = solVector;
    solution u_lsol = A.getSolution(u, sv0);
    for(int ip{0}; ip<m_maxIter; ++ip)
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

        A.initSystem(IGdim);
        // Apply this Picard update before testing convergence or the iteration limit.
        A.assemble(u * grad(u_sol));
        for(index_t Mp=0; Mp<identity_mp.dim(); ++Mp){
            gsMultiPatch<> PsiPsitp_temp;
            timer.restart();
            vsolVector = Poisson.L2ProjectScalar(A.rhs().col(Mp));
            slv_time  += timer.stop();
            v_sol.extract(PsiPsitp_temp);
            Psi.patch(0).coefs().col(Mp) = PsiPsitp_temp.patch(0).coefs().col(0);
        }
        setBoundaryControlPointsAlongNormal(Psi);

        if ( l2Res < tolMAE || ip + 1 == m_maxIter ){
            auto L2MAE        = math::sqrt(ev.integral( pow( 1. - (int_uh_0*abs(rho.val()) + int_uh_1)*jac(PP).det()/CoeffDensity,2)  ));
            auto Ddet         = ev.min(jac(PP).det());
            // ! end Picard loop
            gsInfo  << ". L2_res  : "<<std::scientific<<l2Res
                    << ". L2_err  : "<<std::scientific<<L2MAE
                    << ". min(Jac): "<<std::fixed << std::setprecision(2)<<Ddet
                    << ". max(rho): "<<std::fixed << std::setprecision(2)<<maxrho<<".";
            break;
            } //
    }//END for loop
    Psi.computeTopology();
    // ....
    this->MAmapping = Psi;
    timer.stop();
    gsInfo<<" CPU-time : "<<std::scientific<< slv_time   <<"<>\n";
};

/**
 * \brief Multipatch version of buildMultiPatch(): one adaptive mapping of the unit square per patch.
 *
 * \remark Each Picard sweep updates the 1D edge mappings once, synchronizes shared interfaces,
 * then updates the 2D interiors using the current edges as Dirichlet data.
 * Each coordinate of the mapping solves a Poisson-type problem weighted by the (scaled) density, with Picard iterations.
 * \param density One density patch (on the unit square) per geometry patch.
 */
void gsAdaptiveMultiPatchBuilder::buildMultiPatchMMPDE(const gsMultiPatch<> &density) const
{
    const index_t nP = this->initial_mapping.nPatches();
    GISMO_ENSURE(static_cast<index_t>(density.nPatches()) == nP, "One density patch per geometry patch is required.");
    const bool reuseMapping = !MAmapping.empty();
    GISMO_ENSURE(!reuseMapping || static_cast<index_t>(MAmapping.nPatches()) == nP,
                 "Previous mapping must have one patch per geometry patch.");
    gsInfo<<"<> Multipatch Picard iterations";

    gsStopwatch timer;

    gsExprAssembler<> A(1,1);
    A.options().setSwitch("SameElement",false);
    A.setIntegrationElements(this->m_basis);
    gsExprEvaluator<> ev(A);
    geometryMap G = A.getMap(identity_mp);
    space u       = A.getSpace(this->m_basis);
    gsMatrix<> solVectorx, solVectory, prevx, prevy;
    solution u_solx = A.getSolution(u, solVectorx);
    solution u_soly = A.getSolution(u, solVectory);
    solution u_lsol = A.getSolution(u, prevx);
    solution u_lsoly = A.getSolution(u, prevy);

    // ---- scaling of the density to the range [1, m_IntensityMAE] over all patches ----
    real_t rmax = std::numeric_limits<real_t>::lowest();
    real_t rmin = std::numeric_limits<real_t>::max();
    for (index_t n = 0; n < nP; ++n)
    {
        auto r = A.getCoeff(density.patch(n), G);
        rmax   = std::max<real_t>(rmax, ev.max(abs(r.val())));
        rmin   = std::min<real_t>(rmin, ev.min(abs(r.val())));
    }
    real_t int_uh_0 = 0., int_uh_1 = 1.;
    if (rmax - rmin >= 1e-5 && this->m_IntensityMAE > 1.)
    {
        int_uh_0 = (this->m_IntensityMAE-1.)/(rmax-rmin);
        int_uh_1 = (rmax-this->m_IntensityMAE*rmin)/(rmax-rmin);
    }
    const real_t CoeffDensity = int_uh_0*rmax + int_uh_1 + 5.;

    const gsBasis<>& B = this->m_basis.basis(0);
    const gsMatrix<> bsup = B.support();
    const gsMatrix<> grid = B.anchors();
    const gsMatrix<> idCoefs = B.interpolateData(grid, grid)->coefs();

    // Start edges and interiors from the same mapping.
    gsMultiPatch<> Psi;
    if (reuseMapping)
        Psi = MAmapping;
    else
        for (index_t n = 0; n < nP; ++n)
            Psi.addPatch(B.makeGeometry(idCoefs));
    for (index_t n = 0; n < nP; ++n)
        GISMO_ENSURE(Psi.patch(n).coefsSize() == B.size(),
                     "Previous mapping must use the current moving mesh basis.");

    // Each interface is solved on its first side, then copied to its neighbour.
    std::vector<std::array<bool,4> > isSecond(nP), isInterface(nP), adaptEdge(nP);
    for (auto interface : this->initial_mapping.interfaces())
    {
        isSecond[interface.second().patch][interface.second().index()-1] = true;
        isInterface[interface.first().patch][interface.first().index()-1] = true;
        isInterface[interface.second().patch][interface.second().index()-1] = true;
    }

    std::vector<std::array<gsMatrix<>,4> > edge(nP);
    std::vector<std::array<real_t,4> > edgeA0(nP), edgeA1(nP), edgeCmax(nP);
    for (index_t n = 0; n < nP; ++n)
        for (index_t s = 1; s <= 4; ++s)
        {
            const index_t dir = (s <= 2) ? 1 : 0;
            const gsMatrix<index_t> bnd = B.boundary(boxSide(s));
            const bool fixed = this->m_exactGeo && !isInterface[n][s-1];
            edge[n][s-1].resize(bnd.size(), 1);
            for (index_t i = 0; i < bnd.size(); ++i)
                edge[n][s-1](i) = fixed ? idCoefs(bnd(i), dir)
                                        : Psi.patch(n).coefs()(bnd(i), dir);
            adaptEdge[n][s-1] = !fixed && !isSecond[n][s-1];
            if (!adaptEdge[n][s-1])
                continue;

            // Density is fixed during this call, so edge scaling is computed once.
            EdgeFunction rho(density.patch(n), dir, bsup(1-dir, (s % 2 == 1) ? 0 : 1));
            const gsBasis<>& b1 = B.component(dir);
            gsMatrix<> values;
            rho.eval_into(gsPointGrid<real_t>(bsup(dir,0), bsup(dir,1), 20*b1.size()), values);
            const real_t lo = values.cwiseAbs().minCoeff();
            const real_t hi = values.cwiseAbs().maxCoeff();
            edgeA0[n][s-1] = 0.;
            edgeA1[n][s-1] = 1.;
            if (hi - lo >= 1e-5 && this->m_IntensityMAE > 1.)
            {
                edgeA0[n][s-1] = (this->m_IntensityMAE-1.)/(hi-lo);
                edgeA1[n][s-1] = (hi-this->m_IntensityMAE*lo)/(hi-lo);
            }
            edgeCmax[n][s-1] = edgeA0[n][s-1]*hi + edgeA1[n][s-1] + 5.;
        }

    // The two edge stiffness matrices and solvers are reused in every sweep.
    std::array<gsSparseMatrix<>,2> edgeK, edgeKff;
    std::array<gsSparseSolver<>::CGDiagonal,2> edgeSolver;
    for (index_t dir = 0; dir < 2; ++dir)
    {
        const index_t size = B.component(dir).size();
        edgeK[dir] = assembleStiffness(B.component(dir));
        edgeKff[dir] = edgeK[dir].block(1, 1, size-2, size-2);
        edgeSolver[dir].compute(edgeKff[dir]);
        edgeSolver[dir].setTolerance(1e-14);
    }

    // Perform exactly one 1D Picard update from the current edge coefficients.
    auto updateEdge = [&](const index_t n, const index_t s) -> real_t
    {
        const index_t dir = (s <= 2) ? 1 : 0;
        const gsBasis<>& b1 = B.component(dir);
        const index_t size = b1.size();
        gsMatrix<>& c = edge[n][s-1];
        const real_t a0 = edgeA0[n][s-1], a1 = edgeA1[n][s-1];
        EdgeFunction rho(density.patch(n), dir, bsup(1-dir, (s % 2 == 1) ? 0 : 1));
        gsGeometry<>::uPtr phi = b1.makeGeometry(c);

        gsExprAssembler<> A1(1,1);
        A1.setIntegrationDomain(b1.domain());
        A1.options().setSwitch("SameElement",false);
        space u1 = A1.getSpace(b1);
        solution u_sol = A1.getSolution(u1, c);
        geometryMap PP = A1.getMap(*phi);
        auto frho = A1.getCoeff(rho, PP);
        const real_t Cmax = edgeCmax[n][s-1];
        A1.initSystem();
        A1.assemble(grad(u1) * (1.-(a0*abs(frho.val())+a1)/Cmax) * grad(u_sol).tr());

        gsMatrix<> rhs(size-2, 1);
        for (index_t i = 1; i < size-1; ++i)
            rhs(i-1) = A1.rhs()(i) - edgeK[dir].coeff(i,0)*bsup(dir,0)
                                  - edgeK[dir].coeff(i,size-1)*bsup(dir,1);
        gsMatrix<> next(size, 1);
        next(0) = bsup(dir,0);
        next(size-1) = bsup(dir,1);
        next.block(1, 0, size-2, 1) = edgeSolver[dir].solve(rhs);
        const real_t change = (next-c).norm();
        c = next;
        return change;
    };

    // Keep these geometries alive; boundary conditions refer to their coefficients.
    std::vector<gsGeometry<>::uPtr> gx(nP), gy(nP);
    std::vector<gsBoundaryConditions<> > bc_x(nP), bc_y(nP);
    for (index_t n = 0; n < nP; ++n)
    {
        gx[n] = B.makeGeometry(gsMatrix<>(Psi.patch(n).coefs().col(0)));
        gy[n] = B.makeGeometry(gsMatrix<>(Psi.patch(n).coefs().col(1)));
        bc_x[n].setGeoMap(identity_mp);
        bc_y[n].setGeoMap(identity_mp);
        for (index_t s = 1; s <= 4; ++s)
        {
            bc_x[n].addCondition(0, boxSide(s), condition_type::dirichlet, gx[n].get(), 0, true);
            bc_y[n].addCondition(0, boxSide(s), condition_type::dirichlet, gy[n].get(), 0, true);
        }
    }

    // Precompute the lifting operator -K_free,boundary. New Dirichlet values
    // then require only a matrix-vector product, not a stiffness reassembly.
    u.setup(bc_x[0], dirichlet::interpolation, 0);
    const gsDofMapper mapper = u.mapper();
    const gsSparseMatrix<> stiffness = assembleStiffness(B);
    gsSparseEntries<> entries;
    for (index_t col = 0; col < stiffness.outerSize(); ++col)
        for (gsSparseMatrix<>::InnerIterator it(stiffness, col); it; ++it)
            if (mapper.is_free(it.row()) && mapper.is_boundary(it.col()))
                entries.add(mapper.index(it.row()), mapper.bindex(it.col()), -it.value());
    gsSparseMatrix<> boundaryLifting(mapper.freeSize(), mapper.boundarySize());
    boundaryLifting.setFrom(entries);

    std::vector<gsMatrix<> > xs(nP), ys(nP);
    for (index_t n = 0; n < nP; ++n)
    {
        xs[n].resize(mapper.freeSize(), 1);
        ys[n].resize(mapper.freeSize(), 1);
        for (index_t i = 0; i < B.size(); ++i)
            if (mapper.is_free(i))
            {
                xs[n](mapper.index(i)) = Psi.patch(n).coefs()(i, 0);
                ys[n](mapper.index(i)) = Psi.patch(n).coefs()(i, 1);
            }
    }
    gsMultiPatch<> Psix, Psiy;

    // One interior update, using the edge values from the current sweep.
    auto updateInterior = [&](const index_t n) -> real_t
    {
        geometryMap PP = A.getMap(Psi.patch(n));
        auto frho = A.getCoeff(density.patch(n), PP);

        u.setup(bc_x[n], dirichlet::interpolation, 0);
        A.initSystem();
        solVectorx = xs[n];
        A.assemble((1.-(int_uh_0*abs(frho.val())+int_uh_1)/CoeffDensity) * grad(u) * grad(u_solx).tr());
        prevx = xs[n];
        solVectorx = this->PoissonDir.solve(gsMatrix<>(A.rhs() + boundaryLifting*u.fixedPart()));
        xs[n] = solVectorx;
        real_t res = math::sqrt(ev.integral((u_lsol-u_solx).sqNorm()));
        u_solx.extract(Psix);

        u.setup(bc_y[n], dirichlet::interpolation, 0);
        A.initSystem();
        solVectory = ys[n];
        A.assemble((1.-(int_uh_0*abs(frho.val())+int_uh_1)/CoeffDensity) * grad(u) * grad(u_soly).tr());
        prevy = ys[n];
        solVectory = this->PoissonDir.solve(gsMatrix<>(A.rhs() + boundaryLifting*u.fixedPart()));
        ys[n] = solVectory;
        res = std::max(res, math::sqrt(ev.integral((u_lsoly-u_soly).sqNorm())));
        u_soly.extract(Psiy);

        Psi.patch(n).coefs().col(0) = Psix.patch(0).coefs();
        Psi.patch(n).coefs().col(1) = Psiy.patch(0).coefs();
        return res;
    };

    real_t res = 0., edgeRes = 0.;
    index_t ip = 0;
    while (ip < m_maxIter)
    {
        // 1. Move each independent edge once.
        edgeRes = 0.;
        for (index_t n = 0; n < nP; ++n)
            for (index_t s = 1; s <= 4; ++s)
                if (adaptEdge[n][s-1])
                    edgeRes = std::max(edgeRes, updateEdge(n, s));

        // 2. Give both sides of an interface the same edge distribution.
        for (auto interface : this->initial_mapping.interfaces())
        {
            const patchSide fir = interface.first(), sec = interface.second();
            const index_t dir = (fir.index() <= 2) ? 1 : 0;
            const gsMatrix<> & cf = edge[fir.patch][fir.index()-1];
            if (interface.dirOrientation(fir, dir))
                edge[sec.patch][sec.index()-1] = cf;
            else
                edge[sec.patch][sec.index()-1] = (bsup(dir,0)+bsup(dir,1)) - cf.colwise().reverse().array();
        }

        // 3. Refresh Dirichlet data and the boundary of the current mapping.
        for (index_t n = 0; n < nP; ++n)
            for (index_t s = 1; s <= 4; ++s)
            {
                const index_t dir = (s <= 2) ? 1 : 0;
                const real_t normal = bsup(1-dir, (s % 2 == 1) ? 0 : 1);
                const gsMatrix<index_t> bnd = B.boundary(boxSide(s));
                for (index_t i = 0; i < bnd.size(); ++i)
                {
                    gx[n]->coefs()(bnd(i)) = (dir == 0) ? edge[n][s-1](i) : normal;
                    gy[n]->coefs()(bnd(i)) = (dir == 1) ? edge[n][s-1](i) : normal;
                    Psi.patch(n).coefs()(bnd(i), 0) = gx[n]->coefs()(bnd(i));
                    Psi.patch(n).coefs()(bnd(i), 1) = gy[n]->coefs()(bnd(i));
                }
            }

        // 4. Move each interior once; stop only when edges and interiors settle.
        res = 0.;
        for (index_t n = 0; n < nP; ++n)
            res = std::max(res, updateInterior(n));
        ++ip;
        if (std::max(res, edgeRes) < tolMAE)
            break;
    }

    real_t minJac = std::numeric_limits<real_t>::max();
    for (index_t n = 0; n < nP; ++n)
    {
        geometryMap PP = A.getMap(Psi.patch(n));
        minJac = std::min<real_t>(minJac, ev.min(jac(PP).det()));
    }
    gsInfo << ". Niter: " << ip << ". L2_res: " << std::scientific << res
           << ". edge_res: " << edgeRes
           << ". min(Jac): " << std::fixed << std::setprecision(2) << minJac << ".";

    Psi.interfaces() = this->initial_mapping.interfaces();
    Psi.boundaries() = this->initial_mapping.boundaries();
    this->MAmapping = Psi;
    gsInfo << " CPU-time : " << std::scientific << timer.stop() << "<>\n";
}

/**
 * \brief computes the projection of a composition and return a MultiPatch object :: Collocation
 *
 */
gsMultiPatch<> gsAdaptiveMultiPatchBuilder::buildColCompMultiPatch() const 
{

    gsInfo<<"<Col> computes composition";

    gsMultiPatch<> Psi; 
    double slv_time(0);
    gsStopwatch timer;
    timer.restart();

    for (size_t n = 0; n < this->initial_mapping.nPatches(); ++n)
    {
        const gsMatrix<> initialGrid = mapping_basis.basis(n).anchors();
        // Evaluate f at the Greville points
        gsMatrix<> intervalues         = this->MAmapping.patch(n).eval(initialGrid);
        intervalues                    = intervalues.cwiseMax(0).cwiseMin(1);
        gsMatrix<> finalValues         = this->initial_mapping.patch(n).eval(intervalues);
        gsGeometry<>::uPtr interpolant = mapping_basis.basis(n).interpolateData(finalValues, initialGrid);
        if (this->m_exactGeo)
        {
            // Refine the initial patch until it has as many coefficients as the interpolant, then copy its boundary
            gsGeometry<>::uPtr exact = this->initial_mapping.patch(n).clone();
            GISMO_ENSURE(exact->parDim() == interpolant->parDim(), "Exact geometry and interpolant must have the same parametric dimension.");
            // Preserve the initial degree in every direction; only h-refine below.
            for (index_t d = 0; d < exact->parDim(); ++d)
                GISMO_ENSURE(exact->basis().degree(d) == interpolant->basis().degree(d),
                             "exactGeo: composition must keep the initial degree in direction " << d << ".");
            while (exact->coefsSize() < interpolant->coefsSize())
                exact->uniformRefine();
            if (exact->coefsSize() == interpolant->coefsSize())
            {
                for (gsMultiPatch<>::const_biterator bit = this->initial_mapping.bBegin(); bit != this->initial_mapping.bEnd(); ++bit)
                {
                    if (static_cast<size_t>(bit->patch) != n)
                        continue;
                    const gsMatrix<index_t> bnd = interpolant->basis().boundary(bit->side());
                    for (index_t i = 0; i < bnd.size(); ++i)
                        interpolant->coefs().row(bnd(i)) = exact->coefs().row(bnd(i));
                }
            }
            else
                gsWarn << "exactGeo: refined initial patch " << n << " does not match the basis, boundary not corrected.\n";
        }
        // extract the mapping
        Psi.addPatch(give(interpolant));
    }
    if (this->initial_mapping.nPatches() == 1)
        Psi.computeTopology();
    else
    {
        Psi.interfaces() = this->initial_mapping.interfaces();
        Psi.boundaries() = this->initial_mapping.boundaries();
        Psi.computeTopology();
    }
    //...
    slv_time += timer.stop();
    timer.stop();
    gsInfo<<" CPU-time "<< slv_time   <<"<>\n";
    return Psi;
};

gsMultiPatch<> gsAdaptiveMultiPatchBuilder::buildH1CompMultiPatch(const real_t gradientWeight, const real_t quadValue) const
{
    GISMO_ENSURE(std::isfinite(gradientWeight) && gradientWeight >= 0., "Gradient weight must be finite and non-negative.");
    GISMO_ENSURE(std::isfinite(quadValue) && quadValue > 0., "Quadrature factor must be finite and positive.");

    // Keep the boundary traces, including the existing exactGeo correction.
    gsMultiPatch<> Psi = buildColCompMultiPatch();
    gsBoundaryConditions<> bc;
    for (auto bit = initial_mapping.bBegin(); bit != initial_mapping.bEnd(); ++bit)
        bc.addCondition(*bit, condition_type::dirichlet, nullptr);

    // Each physical coordinate solves (M + alpha K)c = b_value + alpha b_gradient.
    for (index_t c = 0; c < initial_mapping.geoDim(); ++c)
    {
        gsMultiPatch<> target;
        for (size_t p = 0; p < initial_mapping.nPatches(); ++p)
        {
            target.addPatch(gsComposedGeometry<>(MAmapping.patch(p), initial_mapping.patch(p)));
            target.patch(p).coefs() = initial_mapping.patch(p).coefs().col(c);
        }

        gsExprAssembler<> A(1,1);
        A.options().setReal("quA", quadValue);
        A.options().setSwitch("SameElement", false); // Psi may cross knot spans of F.
        A.setIntegrationElements(mapping_basis);
        auto v = A.getSpace(mapping_basis);
        auto f = A.getCoeff(target); // Includes D(F o Psi) through the chain rule.
        v.setup(bc, dirichlet::user, 0);
        const auto & mapper = v.mapper();
        for (size_t p = 0; p < Psi.nPatches(); ++p)
            for (index_t i = 0; i < Psi.basis(p).size(); ++i)
                if (mapper.is_boundary(i, p))
                    v.fixedPart()(mapper.bindex(i, p)) = Psi.patch(p).coefs()(i, c);

        A.initSystem();
        if (A.numDofs() == 0)
            continue;
        A.assemble(v * v.tr() + gradientWeight * grad(v) * grad(v).tr(),
                   v * f + gradientWeight * grad(v) * grad(f).tr());
        gsSparseSolver<>::SimplicialLDLT h1Solver;
        h1Solver.compute(A.matrix());
        GISMO_ENSURE(h1Solver.info() == gsEigen::Success, "H1 composition factorization failed.");
        gsMatrix<> coefficients = h1Solver.solve(A.rhs());
        GISMO_ENSURE(h1Solver.info() == gsEigen::Success && coefficients.allFinite(), "H1 composition solve failed.");
        auto approximation = A.getSolution(v, coefficients);
        for (size_t p = 0; p < Psi.nPatches(); ++p)
        {
            gsMatrix<> values;
            approximation.extract(values, p);
            Psi.patch(p).coefs().col(c) = values.col(0);
        }
    }
    gsInfo << "<H1> composition gradient weight = " << gradientWeight << "\n";
    return Psi;
}


//........................................................................
//........... useful functions for moving meshes B-spline basis ..........
//........................................................................
// Stiffness matrix of a (1D) basis
gsSparseMatrix<> gsAdaptiveMultiPatchBuilder::assembleStiffness(const gsBasis<>& basis) const
{
    gsExprAssembler<> stiff(1,1);
    stiff.setIntegrationDomain(basis.domain());
    typename gsExprAssembler<>::space u = stiff.getSpace(basis);
    stiff.initMatrix();
    stiff.assemble( grad(u) * grad(u).tr() );
    gsSparseMatrix<> result;
    stiff.matrix_into(result);
    return result;
}

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
void gsAdaptiveMultiPatchBuilder::setBoundaryControlPointsAlongNormal(gsMultiPatch<>& Psi) const
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

// Both routines use the same assembly; only the Jacobian dimension differs.
void gsAdaptiveMultiPatchBuilder::assemble_rhsvector_2d(const gsBasis<>& basis,
                           const gsMatrix<double>& vector_mp, const gsMatrix<double>& vector_un,
                           gsMatrix<double>& rhs) const
{
    assembleDensityRhs<2>(basis, vector_mp, vector_un, rhs);
}

void gsAdaptiveMultiPatchBuilder::assemble_rhsvector_3d(const gsBasis<>& basis,
                           const gsMatrix<double>& vector_mp, const gsMatrix<double>& vector_un,
                           gsMatrix<double>& rhs) const
{
    assembleDensityRhs<3>(basis, vector_mp, vector_un, rhs);
}




// // find the knot span index for a given x value
// index_t gsAdaptiveMultiPatchBuilder::find_span(const gsKnotVector<double>& knots, const index_t& degree, const double& x)  const 
// {
//     index_t low = degree;
//     index_t high = knots.size() - 1 - degree;
//     //  Check if point is exactly on left/right boundary, or outside domain
//     if (x <= knots[low]) return low;
//     if (x >= knots[high]) return high - 1;
//     // Binary search for the knot span index
//     index_t mid = (low + high) / 2;
//     while (x < knots[mid] || x >= knots[mid + 1]) {
//         if (x < knots[mid]) high = mid;
//         else low = mid;
//         mid = (low + high) / 2;
//     }
//     return mid;
// };

// // compute the basis functions for a given knot span and degree and x value
// void gsAdaptiveMultiPatchBuilder::basis_functions(const gsKnotVector<double>& knots, const index_t& degree, const double& x, index_t& span,
//                      gsVector<double>& dbasis0,
//                      gsVector<double>& dbasis1)  const 
// {
//     index_t minm = std::min(degree,2);// if degree == 1, then all second derivatives are 0

//     gsVector<double> left(degree), right(degree);
//     gsMatrix<double> ndu(degree + 1, degree + 1);
//     gsMatrix<double> a(2, degree + 1);
//     gsMatrix<double> ders(3, degree + 1);
//     span     = gsAdaptiveMultiPatchBuilder::find_span(knots, degree, x);

//     ndu(0,0)        = 1.0;
//     // computes basis functions
//     // #pragma omp parallel for
//     for (index_t j = 0; j < degree; ++j) {
//         left(j)      = x - knots[span - j];
//         right(j)     = knots[span + j + 1] - x;
//         double saved = 0.0;
//         for (index_t r = 0; r <= j; ++r) {
//             // compute inverse of knot differences and save them into lower triangular part of ndu
//             ndu(j + 1,r)  = 1./(right(r) + left(j - r));
//             //compute basis functions and save them into upper triangular part of ndu
//             double temp   = ndu(r,j) * ndu(j + 1,r);
//             ndu(r,j + 1)  = saved + right(r) * temp;
//             saved         = left(j - r) * temp;
//         }
//         ndu(j + 1, j + 1) = saved;
//     }
//     // compute the derivatives of the basis functions
//     // #pragma omp parallel for
//     for (index_t r = 0; r <= degree; ++r) {
//         ders(0,r) = ndu(r, degree);

//         index_t s1 = 0, s2 = 1;
//         a(0, 0)    = 1.0;
//         for (index_t k = 1; k <= minm; ++k) {
//             double d   = 0.0;
//             index_t rk = r - k;
//             index_t pk = degree - k;
//             if (r >= k) {
//                 a(s2, 0) = a(s1, 0) * ndu(pk + 1, rk);
//                 d        = a(s2, 0) * ndu(rk, pk);
//             }
//             index_t j1 = (rk > -1) ? 1 : -rk;
//             index_t j2 = (r - 1 <= pk) ? k - 1 : degree - r;
//             for (index_t j = j1; j <= j2; ++j)
//                 a(s2,j) = (a(s1, j) - a(s1, j - 1)) * ndu(pk + 1, rk + j);
//             for (index_t j = j1; j <= j2; ++j)
//                 d += a(s2, j) * ndu(rk + j, pk);
//             if (r <= pk) {
//                 a(s2, k) = -a(s1, k - 1) * ndu(pk + 1, r);
//                 d += a(s2, k) * ndu(r, pk);
//             }
//             ders(k, r) = d;
//             std::swap(s1, s2);
//         }
//     }
//     for (index_t j = 0; j <= degree; ++j) {
//         dbasis0(j) = ders(0,j);
//         dbasis1(j) = ders(1,j) * degree;
//     }
// };

// // compute the right-hand side vector for the composition of Left mapping and Monge-Ampere mapping
// void gsAdaptiveMultiPatchBuilder::assemble_rhsvector_1d(
//     const gsBasis<>& basis,
//     const index_t& p1,
//     const gsKnotVector<double>& knots_1,
//     const double& valpt,
//     const index_t& dir_valpt,
//     const index_t& comp_nb,
//     gsVector<double>& rhs) const 
// {
//     rhs.setZero();

//     const double pi = 3.141592653589793;
//     index_t m = 8*p1 + 1;
//     index_t nRoots = (m + 1) / 2;

//     static gsVector<double> u, w;
//     static bool initialized = false;

//     if (!initialized)
//     {
//         u.resize(m);
//         w.resize(m);

//         for (index_t i = 0; i < nRoots; ++i)
//         {
//             double t = std::cos(pi*(i + 0.75)/(m + 0.5));

//             for (index_t j = 0; j < 30; ++j)
//             {
//                 double p0 = 1.0, p1_ = t;

//                 for (index_t k = 1; k < m; ++k)
//                 {
//                     double pn = ((2.0*k+1.0)*t*p1_ - k*p0)/(k+1.0);
//                     p0 = p1_;
//                     p1_ = pn;
//                 }

//                 double denom = (1.0 - t*t);
//                 if (std::abs(denom) < 1e-14) break;

//                 double dp = m*(p0 - t*p1_)/denom;
//                 double dt = -p1_/dp;

//                 t += dt;

//                 if (std::abs(dt) < 1e-14)
//                 {
//                     u(i) = t;
//                     u(m-i-1) = -t;
//                     w(i) = 2.0/((1.0 - t*t)*dp*dp);
//                     w(m-i-1) = w(i);
//                     break;
//                 }
//             }
//         }
//         initialized = true;
//     }

//     index_t nb1 = knots_1.size() - p1 - 1;
//     index_t ne1 = nb1 - p1;

//     gsVector<double> points_cmp(2);
//     points_cmp(dir_valpt) = valpt;
//     index_t dir_op = 1 - dir_valpt;

//     for (index_t ie1 = 0; ie1 < ne1; ++ie1)
//     {
//         double a1 = knots_1[ie1 + p1];
//         double b1 = knots_1[ie1 + p1 + 1];

//         if (std::abs(b1 - a1) < 1e-14)
//             continue;

//         double c0 = 0.5*(a1 + b1);
//         double c1 = 0.5*(b1 - a1);

//         for (index_t g1 = 0; g1 < m; ++g1)
//         {
//             double x1 = c1*u(g1) + c0;
//             double w1 = c1*w(g1);

//             index_t span1;
//             gsVector<double> N(p1+1), dN(p1+1);

//             basis_functions(knots_1, p1, x1, span1, N, dN);

//             points_cmp(dir_op) = x1;
//             // gsInfo << "before :" << points_cmp;
//             gsMatrix<> val_i = MAmapping.patch(0).eval(points_cmp);
//             val_i(dir_valpt) = valpt;
//             val_i.cwiseMax(0).cwiseMin(1);
//             // gsInfo << "after :" << val_i << "---\n";
//             double val_f = initial_mapping.patch(0).eval(val_i)[comp_nb];

//             for (index_t i = 0; i <= p1; ++i)
//             {
//                 index_t gi = span1 - p1 + i;
//                 rhs[gi] += val_f * N(i) * w1;
//             }
//         }
//     }
// };
// // End of the function

// /**
//  * \brief Compute the L^2-projection of a composition and return it as a MultiPatch object.
//  *
//  * \remark The composition is projected onto a B-spline space using an L^2 projection.
//  *         In 2D, a boundary correction can be applied if the approximation near the
//  *         boundary is not sufficiently accurate (currently supported only in the tensor-product case).
//  *
//  * \param quadValue Number of quadrature points used to evaluate the composition.
//  *                  A recommended choice is degree(initial) * degree(MA mapping) + 1.
//  * \param sepBoundary If true, an additional correction is applied at the boundary.
//  */
// gsMultiPatch<> gsAdaptiveMultiPatchBuilder::buildCompMultiPatch(const int quadValue, const bool& sepBoundary) const 
// {
//     if (this->initial_mapping.nPatches() > 1)
//     {
//         gsInfo<<"<Fit> multi-patch case not implemented, calling collocation method instead.\n";
//         return gsAdaptiveMultiPatchBuilder::buildColCompMultiPatch();
//     }
//     gsInfo<<"<L2> computes composition";

//     // target mapping
//     gsMultiPatch<> Psi;
    
//     /// start timing
//     double slv_time(0);
//     gsStopwatch timer;
//     timer.restart();

//     //! [Problem setup] assembler 
//     gsExprAssembler<> A(1,1);
//     A.options().setReal("quA", quadValue);
//     A.options().setSwitch("SameElement",false); // Very important for the composition of the two mappings
//     //...
//     A.setIntegrationElements(mapping_basis);

//     //... 
//     space v        = A.getSpace(mapping_basis);
//     gsMatrix<> vsolVector;
//     solution v_sol = A.getSolution(v, vsolVector);

//     //...
//     geometryMap PP = A.getMap(this->MAmapping);
//     //...
//     A.initSystem();
//     A.assemble(v*v.tr());//Matrix for one component
//     solver.compute( A.matrix() );
//     // ...
//     auto comp      = A.getCoeff(this->initial_mapping, PP);

//     A.initSystem(this->initial_mapping.geoDim());
//     A.assemble(v * comp.tr() );// blocked by this one
//     vsolVector     = solver.solve(A.rhs().col(0));
//     v_sol.extract(Psi);

//     #pragma omp for
//     for(index_t i=1; i<this->initial_mapping.geoDim(); ++i)
//     {
//         gsMultiPatch<> PsiVec;
//         vsolVector = solver.solve(A.rhs().col(i));
//         v_sol.extract(PsiVec);
//         Psi.embed(i+1);
//         Psi.patch(0).coefs().col(i) = PsiVec.patch(0).coefs();
//     }

//     // Correction step
//     if (Psi.dim() == 2 && sepBoundary){
//     // since the initial mapping will be used in the composition.  (F o Psi)
//     const gsKnotVector<double> kv1 =  static_cast<gsTensorNurbs<2> &>( Psi.patch(0)).knots(0);
//     const gsKnotVector<double> kv2 =  static_cast<gsTensorNurbs<2> &>( Psi.patch(0)).knots(1);
//     const index_t degree1 =  static_cast<gsTensorNurbs<2> &>( Psi.patch(0)).degree(0);
//     const index_t degree2 =  static_cast<gsTensorNurbs<2> &>( Psi.patch(0)).degree(1);
//     //...
//     // gsFunctionExpr<> sDir = gsFunctionExpr<>("0",2);
//     // Initialize boundary conditions for Monge-Ampere solver
//     // gsBoundaryConditions<> bc;
//     // bc.setGeoMap(initial_mapping);
//     // for ( gsMultiPatch<>::const_biterator
//     //         bit = initial_mapping.bBegin(); bit != initial_mapping.bEnd(); ++bit)
//     // {
//     //     bc.addCondition( *bit, condition_type::dirichlet, &sDir );
//     // }    
//     gsSparseMatrix<> Mx  = assembleMass(Psi.basis(0).component(0));
//     // eliminateDirichlet1D(boundaryConditionsForDirection(bc,0), A.options(), Mx);
//     gsSparseMatrix<> My  = assembleMass(Psi.basis(0).component(1));
//     // eliminateDirichlet1D(boundaryConditionsForDirection(bc,1), A.options(), My);
//     // componenet by component
//     #pragma omp parallel for
//     for(index_t comp_nb = 0; comp_nb < Psi.geoDim(); ++comp_nb){
//         // y = 0 
//         gsVector<> rhsxy0(Mx.rows());
//         assemble_rhsvector_1d(Psi.basis(0).component(0), degree1, kv1, 0., 1, comp_nb, rhsxy0);
//         solver.compute( Mx );
//         solver.setTolerance(1e-30);
//         auto xsoly0 = solver.solve(rhsxy0);
//         // y = 1
//         gsVector<> rhsxy1(Mx.rows());
//         assemble_rhsvector_1d(Psi.basis(0).component(0), degree1, kv1, 1., 1, comp_nb, rhsxy1);
//         solver.compute( Mx );
//         solver.setTolerance(1e-30);
//         auto xsoly1 = solver.solve(rhsxy1);    
//         // x = 0 
//         gsVector<> rhsx0y(My.rows());
//         assemble_rhsvector_1d(Psi.basis(0).component(1), degree2, kv2, 0., 0, comp_nb, rhsx0y);
//         solver.compute( My );
//         solver.setTolerance(1e-30);
//         auto x0soly = solver.solve(rhsx0y);
//         // x = 1 
//         gsVector<> rhsx1y(My.rows());
//         assemble_rhsvector_1d(Psi.basis(0).component(1), degree2, kv2, 1., 0, comp_nb, rhsx1y);
//         solver.compute( My );
//         solver.setTolerance(1e-30);
//         auto x1soly = solver.solve(rhsx1y);
//         CorrectBoundary(Psi, 0, comp_nb, xsoly0, xsoly1, x0soly, x1soly, true);
//     }
//     gsInfo << ".";
//     }
//     Psi.computeTopology();
//     //...
//     slv_time += timer.stop();
//     timer.stop();
//     gsInfo<<" CPU-time "<< slv_time   <<"<>\n";
//     return Psi;
// };

// gsBoundaryConditions<> gsAdaptiveMultiPatchBuilder::boundaryConditionsForDirection( const gsBoundaryConditions<>& bc, index_t direction ) const
// {
//     gsBoundaryConditions<> result;

//     for ( index_t i = 1; i <= 2; ++i)
//     {
//         patchSide global(0,i+2*direction), local(0,i);
//         const boundary_condition<double>* cond = bc.getConditionFromSide(global);
//         if (cond)
//             result.addCondition(local,cond->type(),cond->function());
//     }
//     return result;
// }


// // Correct the control point at the boundary of a final mapping 
// void gsAdaptiveMultiPatchBuilder::CorrectBoundary(gsMultiPatch<>& Psi, const index_t& patchNumber, const index_t& patch_cmp, const gsMatrix<>& xsoly0, const gsMatrix<>& xsoly1, const gsMatrix<>& x0soly, const gsMatrix<>& x1soly, const bool& corners) const
// {
//     // ...
//     if (corners){
//     // update corners
//     for (int i_x =0; i_x < Psi.patch(patchNumber).basis().boundary(1).size(); ++i_x)
//     {
//     // x=0 control points be like (0,:) in this case
//     Psi.patch(patchNumber).coef( Psi.patch(patchNumber).basis().boundary(1).at(i_x) ).array()[patch_cmp] = x0soly(i_x);
//     // x=1 control points be like (1,:) in this case
//     Psi.patch(patchNumber).coef( Psi.patch(patchNumber).basis().boundary(2).at(i_x) ).array()[patch_cmp] = x1soly(i_x);
//     }

//     for (int i_x =0; i_x < Psi.patch(patchNumber).basis().boundary(3).size(); ++i_x) 
//     {
//     // y=0 control points be like (:,0) in this case
//     Psi.patch(patchNumber).coef( Psi.patch(patchNumber).basis().boundary(3).at(i_x) ).array()[patch_cmp] = xsoly0(i_x);
//     // y=1 control points be like (:,1) in this case
//     Psi.patch(patchNumber).coef( Psi.patch(patchNumber).basis().boundary(4).at(i_x) ).array()[patch_cmp] = xsoly1(i_x);
//     }
//     }
//     else{
//     for (int i_x =0; i_x < Psi.patch(patchNumber).basis().boundary(1).size()-1; ++i_x)
//     {
//     // x=0 control points be like (0,:) in this case
//     Psi.patch(patchNumber).coef( Psi.patch(patchNumber).basis().boundary(1).at(i_x+1) ).array()[patch_cmp] = x0soly(i_x);
//     // x=1 control points be like (1,:) in this case
//     Psi.patch(patchNumber).coef( Psi.patch(patchNumber).basis().boundary(2).at(i_x+1) ).array()[patch_cmp] = x1soly(i_x);
//     }

//     for (int i_x =0; i_x < Psi.patch(patchNumber).basis().boundary(3).size()-1; ++i_x) 
//     {
//     // y=0 control points be like (:,0) in this case
//     Psi.patch(patchNumber).coef( Psi.patch(patchNumber).basis().boundary(3).at(i_x+1) ).array()[patch_cmp] = xsoly0(i_x);
//     // y=1 control points be like (:,1) in this case
//     Psi.patch(patchNumber).coef( Psi.patch(patchNumber).basis().boundary(4).at(i_x+1) ).array()[patch_cmp] = xsoly1(i_x);
//     }        
//     }
// }

// /**
//  * \brief Compute the projection (fitting) of a composition and return it as a MultiPatch object.
//  *
//  * \remark The composition is approximated in a B-spline space using a least-squares (L^2) fitting approach.
//  *         In 2D, a boundary correction can be applied if the approximation near the boundary
//  *         is not sufficiently accurate (currently supported only in the tensor-product case).
//  *
//  * \param numElData Refinement level of the sampling grid used for the least-squares fitting,
//  *                  relative to the identity mapping.
//  * \param lambda Regularization (preconditioning) parameter; recommended value is 10^{-6},
//  *               or at least smaller than 10^{-3}.
//  * \param sepboundary if one want to separate boundary from inner  points in fitting procedure (based on Newton iteration( solve optimization pr))
//  */
// gsMultiPatch<> gsAdaptiveMultiPatchBuilder::buildFitCompMultiPatch(const int numElData, const real_t lambda, const bool& sepboundary) const 
// {
//     if (this->initial_mapping.nPatches() > 1 || mapping_basis.basis(0).weights().any())
//     {
//         gsInfo<<"<Fit> multi-patch case not implemented, calling collocation method instead.\n";
//         return gsAdaptiveMultiPatchBuilder::buildColCompMultiPatch();
//     }
//     gsInfo<<"<Fit> computes composition";
//     assert(mapping_basis.dim() == 2 && "Only single patch 2D fitting is implemented so far.");

//     gsMultiPatch<> Psi;    
//     // Copy tensor basis
//     gsTHBSplineBasis<2>  THB ( mapping_basis.basis(0));

//     double slv_time(0);
//     gsStopwatch timer;
//     timer.restart();
//     //...  just to generate grid in the computational domain  (we start from identity as we want tensor grid)
//     gsMultiBasis<> T_tbasis(identity_mp, true);

//     while ( T_tbasis.basis(0).size() <= std::max(mapping_basis.basis(0).size()*30, identity_mp.basis(0).size()*numElData))
//     {
//         T_tbasis.uniformRefine();
//     }
//     gsInfo<<":gridsize="<<T_tbasis.basis(0).size()<<"/"<<mapping_basis.basis(0).size();
//     gsMatrix<> initialGrid    = T_tbasis.basis(0).anchors();

//     // ------------------------------
//     // ... generate data for fitting
//     // gsMatrix<> initialGrid             = mapping_basis.basis(0).anchors();
//     // Evaluate f at the Greville points
//     gsMatrix<> intervalues          = this->MAmapping.patch(0).eval(initialGrid);
//     intervalues                     = intervalues.cwiseMax(0).cwiseMin(1);
//     gsMatrix<> finalValues           = this->initial_mapping.patch(0).eval(intervalues);


//     //-----------------------
//     // start fitting
//     //-----------------------

//     // Create hierarchical refinement object
//     gsFitting<> ref(initialGrid, finalValues, THB);
//     //... compute coefs
//     ref.compute(lambda);

//     if(sepboundary){
//     // number of grid points in one direction, assuming a square grid.
//     index_t ngrids        = sqrt(initialGrid.cols()); 
//     std::vector<index_t> boundaryIdx;
//     for (index_t i = 0; i < ngrids; ++i)
//     {
//         // y = 0 (bottom)
//         boundaryIdx.push_back(i * ngrids + 0);
//         // y = 1 (top)
//         boundaryIdx.push_back(i * ngrids + (ngrids - 1));
//         // x = 0 (left)
//         boundaryIdx.push_back(0 * ngrids + i);
//         // x = 1 (right)
//         boundaryIdx.push_back((ngrids - 1) * ngrids + i);
//     }
//     // Remove duplicates from boundary indices
//     std::sort(boundaryIdx.begin(), boundaryIdx.end());
//     boundaryIdx.erase(std::unique(boundaryIdx.begin(), boundaryIdx.end()),
//                     boundaryIdx.end());
//     // Number of interior points
//     index_t totalPts = ngrids * ngrids;
//     index_t nbInterior = totalPts - boundaryIdx.size();
//     // Combine interior and boundary indices for fitting
//     std::vector<index_t> interpIdx;
//     interpIdx.push_back(nbInterior);
//     for (auto b : boundaryIdx)
//         interpIdx.push_back(b);
//     //... initialization of geometry for fitting
//     // gsMultiPatch<> PsiCol = initial_mapping;
//     // while( PsiCol.basis(0).numElements() < mapping_basis.basis(0).numElements())   
//     //     PsiCol.uniformRefine();
//     // ref.initializeGeometry(PsiCol.patch(0).coefs(), initialGrid);
//     ref.parameterProjectionSepBoundary(1e-4, interpIdx);
//     }

//     //! [extract the mapping]
//     Psi.addPatch(give(*ref.result()));
//     Psi.computeTopology();
//     //...
//     slv_time += timer.stop();
//     timer.stop();
//     gsInfo<<" CPU-time "<< slv_time   <<"<>\n";
//     return Psi;
// };