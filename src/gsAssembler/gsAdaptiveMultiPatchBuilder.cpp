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

// Constructor implementation
gsAdaptiveMultiPatchBuilder::gsAdaptiveMultiPatchBuilder(const gsMultiBasis<> basis,
                            const gsMultiPatch<> mapping,
                            index_t maxIter,
                            double IntensityMAE,
                            index_t numReduce)
{
    gsInfo<<"\n <>r-refinement (!!!";
    gsMultiBasis<double> dbasis(mapping, true);//not NURBS
    if (dbasis.degree()-numReduce >= 1)
        dbasis.degreeReduce(numReduce);
    //... condition for the convergence 
    while (dbasis.basis(0).numElements()<basis.basis(0).numElements())
        dbasis.uniformRefine();
    gsInfo<<" We use B-spline of degree "<< dbasis.degree()<<" for AdMapping) \n";
    this->m_basis        = dbasis;
    this->n_basis        = basis;
    this->m_mapping      = mapping; 
    this->m_maxIter      = maxIter;
    this->m_IntensityMAE = IntensityMAE;

    // .... one single patch
    auto corners         = basis.basis(0).support();
    mp                   = gsNurbsCreator<>::BSplineSquareGrid(1,1,corners.at(2), corners.at(0), corners.at(1));
    gsFunctionExpr<> sN("x","y",2);
    if (basis.dim() == 3){
        mp  = gsNurbsCreator<>::BSplineCubeGrid(1,1,1,corners.at(3),corners.at(0)-0.5,corners.at(1)-0.5,corners.at(2)-0.5);
        // Manufactured identity mapping
         sN = gsFunctionExpr<>("x","y","z",3);
    }
    //Get all interfaces and boundaries:
    mp.computeTopology();

    // ... Neumann boundary conditions
    gsBoundaryConditions<> bc_mae;
    bc_mae.setGeoMap(mp);
    // For simplicity, set Neumann boundary conditions
    for ( gsMultiPatch<>::const_biterator
            bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
    {
        bc_mae.addCondition( *bit, condition_type::neumann, &sN );
    }
    
    //! [Problem setup]
    gsExprAssembler<> A(1,1);

    //::::::::::::::::::::      Poisson fast diagonalization solver         :::::::::::::::::::::::::
    gsPatchPreconditionersCreator<double>::Poisson_FastDiag Poisson(this->m_basis.basis(0), bc_mae, A.options(), 1e-6);  
    this->Poisson   = Poisson;
    this->mp        = mp;
    gsInfo<<"<> \n";
}

// Method to Project normal control points
void gsAdaptiveMultiPatchBuilder::ProjectionNormalCPoints(gsMultiPatch<>& Psi, int boxMaxNumber) const
{
    // Projection normal of control points (exact geometry)
    for (int boxNumber = 0; boxNumber < boxMaxNumber; ++boxNumber)
    {
        index_t bndIter  = 1;
        #pragma omp parallel for
        for (index_t dim = 0; dim < Psi.dim(); ++dim){
            float_t lVal    = this->mp.patch(boxNumber).coef( this->mp.patch(boxNumber).basis().boundary(bndIter).at(0) ).array()[dim];
            for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(bndIter).size(); ++i_x) // x=0 control points be like (0,:) in this case
            {
                Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(bndIter).at(i_x) ).array()[dim] = lVal;
            }
            bndIter ++;
            float_t hVal    = this->mp.patch(boxNumber).coef( this->mp.patch(boxNumber).basis().boundary(bndIter).at(0) ).array()[dim];
            for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(bndIter).size(); ++i_x)// x=1 control points be like (1,:) in this case
            {
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(bndIter).at(i_x) ).array()[dim] = hVal;
            }
            bndIter ++;
        }
    }
}

// Build and return a density as a MultiPatch object from analytical function
gsMultiPatch<> gsAdaptiveMultiPatchBuilder::buildAnalyticDensity(const gsFunctionExpr<>   &f) const 
{
    gsInfo<<"<>density function";
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    // Elements used for numerical integration
    A.setIntegrationElements(this->m_basis);
    gsExprEvaluator<> ev(A);

    // Set the discretization space
    space u             = A.getSpace(this->m_basis);
    // Set the Target geometry map
    geometryMap GLeft   = A.getMap(this->m_mapping);
    // Set the source term with respect to target geometry
    auto ff             = A.getCoeff(f, GLeft);
    // Solution vector and solution variable
    gsMatrix<> densityVector;

    //u.setup(bc_mae, dirichlet::l2Projection, 0);
    A.initSystem();
    A.assemble(u* ff.val()); //rhs vector
    densityVector        = this->Poisson.L2ProjectScalar(A.rhs());
    //...
    solution density_sol = A.getSolution(u, densityVector);
    gsMultiPatch<> density;
    density_sol.extract(density);
    gsInfo<<"<>\n";
    return  density;
}

// Build and return a density as a MultiPatch object from solution vector
gsMultiPatch<> gsAdaptiveMultiPatchBuilder::buildDensity(const std::vector<double> &elwiseERROR, const double eps, bool maxminVar) const 
{
    gsInfo<<"<>density function";
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;
    //...............error as a piecewise constant function
    gsMultiBasis<double> basis_0(m_mapping, true);
    basis_0.setDegree(0);
    //... refine uniformly until the number of elements is equal to the number of elwiseERROR
    int refnumb = basis_0.basis(0).size();
    int elwnumb = elwiseERROR.size();
    while ( refnumb< elwnumb)
    {
        basis_0.uniformRefine();
        refnumb = basis_0.basis(0).size();
    }
    gsExprAssembler<> A_0(1,1);
    // Elements used for numerical integration
    A_0.setIntegrationElements(basis_0);    
    // Set the discretization space
    space u_0          = A_0.getSpace(basis_0);
    A_0.initSystem();
    auto errorVector   = A_0.rhs();
    solution error_sol = A_0.getSolution(u_0, errorVector);
    // ...
    #pragma omp parallel for
    for (index_t i = 0; i < elwnumb; i++){
        errorVector(i) = elwiseERROR[i];
    }
    //... normalize the error vector
    if (maxminVar){
        const double Maxvalue   = errorVector.maxCoeff();
        const double Minvalue   = errorVector.minCoeff();
        // gsInfo << "Density function: min "<< errorVector.minCoeff() <<"/ max " << errorVector.maxCoeff() << "\n";    
        const double meanvalue  = eps*(Maxvalue + Minvalue);
        for (index_t i1 = 0; i1 < elwnumb; i1++){
            if (errorVector(i1) > Minvalue+meanvalue)
                errorVector(i1)  = Minvalue+meanvalue;
            //else errorVector(i1) = Minvalue;
        }
    }

    //...............End error as a function
    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    // Elements used for numerical integration
    A.setIntegrationElements(this->m_basis);
    gsExprEvaluator<> ev(A);

    // Set the discretization space
    space u             = A.getSpace(this->m_basis);
    // Set the Target geometry map
    geometryMap GLeft   = A.getMap(this->m_mapping);
    // Set the source term with respect to target geometry
    //auto            ff  = A.getCoeff(Multipatchsol, GLeft);
    // Solution vector and solution variable
    gsMatrix<> densityVector;
    //...
    solution density_sol = A.getSolution(u, densityVector);
    //...
    //u.setup(bc_mae, dirichlet::l2Projection, 0);
    A.initSystem();
    A.assemble(u * error_sol); //rhs vector

    densityVector        = this->Poisson.L2ProjectScalar(A.rhs());
    //...
    gsMultiPatch<> density;
    density_sol.extract(density);
    gsInfo<<"<>\n";
    return  density;
}

// Build and return a density as a MultiPatch object from solution vector using local h-refinement strategies
gsMultiPatch<> gsAdaptiveMultiPatchBuilder::buildStrategyDensity(const std::vector<double> &elwiseERROR, const double MarkPercentage) const 
{
    gsInfo<<"<>density function";
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;
    //...............error as a piecewise constant function
    gsMultiBasis<double> basis_0(m_mapping, true);
    basis_0.setDegree(0);
    //... refine uniformly until the number of elements is equal to the number of elwiseERROR
    int refnumb = basis_0.basis(0).size();
    int elwnumb = elwiseERROR.size();
    while ( refnumb< elwnumb)
    {
        basis_0.uniformRefine();
        refnumb = basis_0.basis(0).size();
    }
    gsExprAssembler<> A_0(1,1);
    // Elements used for numerical integration
    A_0.setIntegrationElements(basis_0);    
    // Set the discretization space
    space u_0          = A_0.getSpace(basis_0);
    A_0.initSystem();
    auto errorVector   = A_0.rhs();
    solution error_sol = A_0.getSolution(u_0, errorVector);
    //..
    // The vector of local errors will need to be sorted,
    // which will be done on a copy:
    std::vector<double> elErrCopy = elwiseERROR;

    // Compute the index from which the refinement should start,
    // once the vector is sorted.
    size_t idxRefineStart = cast<double,size_t>( math::floor( MarkPercentage * (double)(elwnumb) ) );
    // ...and just to be sure we are in range:
    if( idxRefineStart == elErrCopy.size() )
    {
        GISMO_ASSERT(idxRefineStart >= 1, "idxRefineStart can't get negative");
        idxRefineStart -= 1;
    }
    // Sort the list using bubblesort.
    // After each loop, the largest elements are at the end
    // of the list. Since we are only interested in the largest elements,
    // it is enough to run the sorting until enough "largest" elements
    // have been found, i.e., until we have reached indexRefineStart
    size_t lastSwapDone = elErrCopy.size() - 1;
    size_t lastCheckIdx = lastSwapDone;

    bool didSwap;
    double tmp;
    do{
        didSwap = false;
        lastCheckIdx = lastSwapDone;
        for( size_t i=0; i < lastCheckIdx; i++)
            if( elErrCopy[i] > elErrCopy[i+1] )
            {
                tmp = elErrCopy[i];
                elErrCopy[i] = elErrCopy[i+1];
                elErrCopy[i+1] = tmp;

                didSwap = true;
                lastSwapDone = i;
            }
    }while( didSwap && (lastSwapDone+1 >= idxRefineStart ) );

    // Compute the threshold:
    auto Thr = elErrCopy[ idxRefineStart ];
    // Now just check for each element, whether the local error
    // is above the computed threshold or not, and mark accordingly.
    for( size_t i=0; i < elwiseERROR.size(); i++)
        ( elwiseERROR[i] >= Thr ? errorVector(i) = Thr : errorVector(i) = elwiseERROR[i] );

    //...............End error as a function
    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    // Elements used for numerical integration
    A.setIntegrationElements(this->m_basis);
    gsExprEvaluator<> ev(A);

    // Set the discretization space
    space u             = A.getSpace(this->m_basis);
    // Solution vector and solution variable
    gsMatrix<> densityVector;
    //...
    solution density_sol = A.getSolution(u, densityVector);
    //...
    //u.setup(bc_mae, dirichlet::l2Projection, 0);
    A.initSystem();
    A.assemble(u * error_sol); //rhs vector

    densityVector        = this->Poisson.L2ProjectScalar(A.rhs());
    //...
    gsMultiPatch<> density;
    density_sol.extract(density);
    gsInfo<<"<>\n";
    return  density;
}


// Build and return a MultiPatch object
gsMultiPatch<> gsAdaptiveMultiPatchBuilder::buildMultiPatch(const gsMultiPatch<> &density) const 
{
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Neumann Boundary conditions object to define and manage boundary conditions for the problem
    gsBoundaryConditions<> bc_mae;
    // Manufactured identity mapping
    gsFunctionExpr<> sN("x","y",2);
    if (this->m_basis.dim() == 3){
        sN = gsFunctionExpr<>("x","y","z",3);
    }

    bc_mae.setGeoMap(this->mp);
    // For simplicity, set Neumann boundary conditions
    for ( gsMultiPatch<>::const_biterator
                bit = this->mp.bBegin(); bit != this->mp.bEnd(); ++bit)
    {
            bc_mae.addCondition( *bit, condition_type::neumann, &sN );
    }
    
    // Target mapping
    gsMultiPatch<> Psi;
    gsInfo<<"<> Picard iterations \n";
    double slv_time(0);
    gsStopwatch timer;
    timer.restart();
    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    // It could be beneficial for the composition of the two mappings
    // A.options().setReal("quA", quadValue);
    //A.options().setInt("quB", 2);
    A.options().setSwitch("SameElement",false);

    // Elements used for numerical integration
    A.setIntegrationElements(this->m_basis);
    gsExprEvaluator<> ev(A);

    // Set the geometry map
    geometryMap G  = A.getMap(this->mp);
    // Set the discretization space
    space u        = A.getSpace(this->m_basis);

    // Set pow for BFO method dim in parameteric domain
    auto IGdim     = G.domainDim();

    // Set factor for BFO method
    auto gammaMAE  = factorial(G.domainDim());

    // Solution vector and solution variable
    gsMatrix<> solVector;
    solution u_sol   = A.getSolution(u, solVector);

    // Solution vector and solution for projecting each direction of mapping
    gsMatrix<> vsolVector;
    solution v_sol = A.getSolution(u, vsolVector);

    // ---- manipulation of density function ----
    auto rho         = A.getCoeff(density);
    auto empldensity = (ev.max(abs(rho.val()))-ev.min(abs(rho.val())));
    double  int_uh_0 = 0.;
    double  int_uh_1 = 1.;
    if (empldensity < 1e-5|| this->m_IntensityMAE <= 1. )
    {
        gsInfo << " rho = 1.\n";
    }
    else{
        int_uh_0  = (this->m_IntensityMAE-1.)/empldensity;
        int_uh_1  = (1.*ev.max(abs(rho.val()))-this->m_IntensityMAE*ev.min(abs(rho.val())))/empldensity;
        //gsInfo << "Density function is not constant in the domain\n";
    gsInfo << "rho = "<< ev.min(int_uh_0*abs(rho.val()) + int_uh_1)<<"/" << ev.max(int_uh_0*abs(rho.val()) + int_uh_1) << std::flush;
    }
    // ......... End initialization for density.........
    //u.setup(bc_mae, dirichlet::l2Projection, 0);
    // Compute the system matrix and right-hand side

    // Initialize the system :  identity mapping as initial guess
    A.initSystem();
    // Initialize the system : start Computing the conductivity coeffeicient ...
    // Compute the Neumann terms defined on physical space
    auto g_N      = A.getBdrFunction(G);
    auto Neumann_Int{ev.integralBdrBc(bc_mae.get("Neumann"), g_N.tr() * nv(G) )};
    //... nromalisation of density function
    auto CoeffDensity{ev.integral((int_uh_0*abs(rho.val()) + int_uh_1))};
    auto ExprMAE     = pow(pow(IGdim,IGdim)-gammaMAE+gammaMAE * CoeffDensity/(int_uh_0*abs(rho.val()) + int_uh_1), 1./IGdim);
    auto CoeffConductivity{Neumann_Int/ev.integral(ExprMAE)};

    A.assemble(u*  CoeffConductivity * (-1.)*ExprMAE); //rhs vector

    // Compute the Neumann terms defined on physical space
    A.assembleBdr(bc_mae.get("Neumann"), u * g_N.tr() * nv(G) );
    //A.assembleIfc(mp.interfaces(), u.left() * (u_I.tr() * nv(G.left())));
    //A.assembleIfc(mp.interfaces(), u.right() * (u_I.tr() * nv(G.right())));

    gsInfo<< "." <<std::flush;// Assemblying done
    solVector = A.rhs();
    solVector = this->Poisson.solve(A.rhs());
    gsInfo<< "." <<std::flush; // Linear solving done

    // Initial guess for the gradient of potential function
    A.initSystem(IGdim);
    // Obtain control points for the gradient of Psi
    A.assemble(u * grad(u_sol));
    timer.restart();
    vsolVector    = Poisson.L2ProjectScalar(A.rhs().col(0));
    slv_time     += timer.stop();
    v_sol.extract(Psi);
    for(index_t Mp=0; Mp<mp.dim(); ++Mp){
    gsMultiPatch<> PsiPsitp_temp;
    timer.restart();
    vsolVector    = Poisson.L2ProjectScalar(A.rhs().col(Mp));
    slv_time     += timer.stop();
    v_sol.extract(PsiPsitp_temp);
    Psi.embed(Mp+1);
    Psi.patch(0).coefs().col(Mp) = PsiPsitp_temp.patch(0).coefs().col(0);
    }
    // // ... correct boundary
    ProjectionNormalCPoints(Psi);
    //! [Solver loop]
    gsInfo<< A.numDofs() <<std::flush;
    // Picard loop
    auto  sv0 = solVector; //
    solution u_lsol = A.getSolution(u, sv0);
    for(int ip{0}; ip<=m_maxIter; ++ip)
    {
        //.. geometry map
        geometryMap PP    = A.getMap(Psi);
        //... density in new optimized mesh
        auto rho = A.getCoeff(density, PP);
        // ... update residual
        solution u_sol    = A.getSolution(u, solVector);

        // ...  0  dirichlet for boundaries
        //u.setup(bc_mae, dirichlet::l2Projection, 0);
        // Initialize the system
        A.initSystem();

        // Compute the system matrix and right-hand side ... Monge-Ampere eqaution .....
        // .. update Coeffeicient of conductivity
        auto  ExprMAE     = pow( abs(pow(div(PP).val(),IGdim) - gammaMAE*jac(PP).det())+ gammaMAE*CoeffDensity/(int_uh_0*abs(rho.val()) + int_uh_1), 1./IGdim);
        // if (dbasis.minCwiseDegree() > 2)
        // auto  ExprMAE     = pow( abs(pow(lapl(u_sol).val(),IGdim) - gammaMAE*hess(u_sol).det())+ gammaMAE*CoeffDensity/(int_uh_0*abs(rho.val()) + int_uh_1), 1./IGdim);
        auto IntegDensity = ev.integral(ExprMAE);
        CoeffConductivity = Neumann_Int/IntegDensity;
        // MAE system
        A.assemble(u * CoeffConductivity * (-1.) * ExprMAE);//rhs vector
        //gsInfo << "End Assemnles \n";

        // Compute the Neumann terms defined on physical space
        auto g_N = A.getBdrFunction(G);
        A.assembleBdr(bc_mae.get("Neumann"), u * g_N.tr() * nv(G) );
        //A.assembleIfc(mp.interfaces(), u.left() * (u_I.tr() * nv(G.left())));
        //A.assembleIfc(mp.interfaces(), u.right() * (u_I.tr() * nv(G.right())));

        gsInfo<< " ." <<std::flush;// Assemblying done
        solVector = this->Poisson.solve(A.rhs());
        gsInfo<< "." <<std::flush; // Linear solving done

        auto l2errRes = math::sqrt(ev.integral( ( grad(u_sol)-grad(u_lsol)).sqNorm()  ));
        auto L2MAERes = math::sqrt(ev.integral( pow( 1. - (int_uh_0*abs(rho.val()) + int_uh_1)*jac(PP).det()/CoeffDensity,2)  ));
        auto Ddet     = ev.min(jac(PP).det());
        sv0           = solVector;

        if ( l2errRes < 1e-8 || ip == m_maxIter ){
            // ! end Picard loop
            gsInfo<< "\n Niter in Picard : " << ip
                    << ".. H1 residual : "<<std::scientific<<l2errRes
                    << ".. L2 MAE residual : "<<std::scientific<<L2MAERes
                    << ".. min JAcobian : "<<Ddet<<"..";
            break;
            } //
    A.initSystem(IGdim);
    // Obtain control points for the gradient of Psi
    A.assemble(u * grad(u_sol));
    for(index_t Mp=0; Mp<mp.dim(); ++Mp){
    gsMultiPatch<> PsiPsitp_temp;
    timer.restart();
    vsolVector = Poisson.L2ProjectScalar(A.rhs().col(Mp));
    slv_time     += timer.stop();
    v_sol.extract(PsiPsitp_temp);
    Psi.patch(0).coefs().col(Mp) = PsiPsitp_temp.patch(0).coefs().col(0);
    }
    //... correct the boundary
    ProjectionNormalCPoints(Psi);
    // ...
    }//for loop

    Psi.addAutoBoundaries();
    Psi.computeTopology();
    timer.stop();
    gsInfo<<" CPU-time : "<< slv_time   <<"<>\n";
    return Psi;

};

// Build and return a MultiPatch object in hierarchical refinement
gsMultiPatch<> gsAdaptiveMultiPatchBuilder::buildHBMultiPatch(const gsMultiBasis<> dbasis, const gsMultiPatch<> &density) const 
{
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Neumann Boundary conditions object to define and manage boundary conditions for the problem
    gsBoundaryConditions<> bc_mae;
    // Manufactured identity mapping
    gsFunctionExpr<> sN("x","y",2);
    if (this->m_basis.dim() == 3){
        sN = gsFunctionExpr<>("x","y","z",3);
    }

    bc_mae.setGeoMap(this->mp);
    // For simplicity, set Neumann boundary conditions
    for ( gsMultiPatch<>::const_biterator
                bit = this->mp.bBegin(); bit != this->mp.bEnd(); ++bit)
    {
            bc_mae.addCondition( *bit, condition_type::neumann, &sN );
    }

    gsInfo<<"<> Picard iterations \n";
    double slv_time(0);
    gsStopwatch timer;
    timer.restart();

    gsMultiPatch<> Psi;

    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    // It could be beneficial for the composition of the two mappings
    // A.options().setReal("quA", quadValue);
    //A.options().setInt("quB", 2);
    A.options().setSwitch("SameElement",false);

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);
    // Set the square geometry map
    geometryMap G     = A.getMap(mp);

    // Set the Target geometry adaptive map
    geometryMap PP    = A.getMap(Psi);

    // Set pow for BFO method dim in parameteric domain
    auto IGdim        = G.domainDim();

    // Set factor for BFO method
    auto gammaMAE     = factorial(G.domainDim());

    // Set the discretization space
    space u           = A.getSpace(dbasis);

    //... solution vector and solution variable for gradient of potential function
    gsMatrix<> vsolVector;
    solution v_sol    = A.getSolution(u, vsolVector);

    //gsFunctionExpr<> sI("0.5*(x**2+y**2)+x*y",2);
    auto u_I          = ev.getVariable(sN, G);

    // Solution vector and solution variable
    gsMatrix<> solVector;
    solution u_sol    = A.getSolution(u, solVector);

    //::::::::::::::::::::      mesh adaptation solver         :::::::::::::::::::::::::
    gsSparseSolver<>::CGDiagonal solver; // exact solver
    gsSparseSolver<>::CGDiagonal Ssolver;// relaxation solver
    // Ssolver.setMaxIterations(20);
    // Ssolver.setTolerance(1e-8);

    timer.restart();
    // -------------------- for projection --------------------
    u.setup(bc_mae, dirichlet::l2Projection, 0);
    A.initSystem();
    A.assemble(u *u.tr());//matrix
    auto MProj = A.matrix();
    solver.compute( MProj);
    // -------------------- for MAE system ----------------
    u.setup(bc_mae, dirichlet::l2Projection, 0);
    A.initSystem();
    A.assemble(grad(u) * grad(u).tr() + 1e-6 * u *u.tr());//matrix
    auto MMAe = A.matrix();
    Ssolver.compute( MMAe );
    slv_time += timer.stop();

    //---------------- renormalize the density function -------------------
    // Solution vector and solution variable
    auto rho_in      = A.getCoeff(density, G);
    auto rho         = A.getCoeff(density, PP);

    // ... manipulation of density function
    auto empldensity = (ev.max(abs(rho_in.val()))-ev.min(abs(rho_in.val())));
    double  int_uh_0 = 0.;
    double  int_uh_1  = 1.;
    if (empldensity < 1e-5|| this->m_IntensityMAE <= 1. )
    {
        gsInfo << "Density function is constant in the domain rho = 1.\n";
    }
    else{
        int_uh_0  = (this->m_IntensityMAE-1.)/empldensity;
        int_uh_1  = (1.*ev.max(abs(rho_in.val()))-this->m_IntensityMAE*ev.min(abs(rho_in.val())))/empldensity;
        gsInfo << "Density functio: min "<< ev.min(int_uh_0*abs(rho_in.val()) + int_uh_1)<<"/ max " << ev.max(int_uh_0*abs(rho_in.val()) + int_uh_1) << "\n";
    }
    // ......... End initialization for density.........
        
    // ......... Start solving the Monge-Ampere equation .........
    A.initSystem();
    // Initialize the system : start Computing the conductivity coeffeicient ...
    // Compute the Neumann terms defined on physical space
    auto g_N      = A.getBdrFunction(G);
    auto Neumann_Int{ev.integralBdrBc(bc_mae.get("Neumann"), g_N.tr() * nv(G) )};
    //... normalisation of a density function
    auto CoeffDensity{ev.integral((int_uh_0*abs(rho_in.val()) + int_uh_1))};
    auto ExprMAE     = pow(pow(IGdim,IGdim)-gammaMAE+gammaMAE * CoeffDensity/(int_uh_0*abs(rho_in.val()) + int_uh_1), 1./IGdim);
    auto CoeffConductivity{Neumann_Int/ev.integral(ExprMAE)};

    A.assemble(u*  CoeffConductivity * (-1.)*ExprMAE//rhs vector
    );
    // Compute the Neumann terms defined on physical space
    A.assembleBdr(bc_mae.get("Neumann"), u * g_N.tr() * nv(G) );
    gsInfo<< "." <<std::flush;// Assemblying done
    timer.restart();
    solVector = A.rhs();
    // solVector = Poisson.solve(A.rhs());
    solVector = Ssolver.solve(A.rhs());
    slv_time += timer.stop();
    gsInfo<< "." <<std::flush; // Linear solving done

    // Initial guess for the gradient of potential function
    A.initSystem(IGdim);
    // Obtain control points for the gradient of Psi
    A.assemble(u * grad(u_sol));
    timer.restart();
    vsolVector    = solver.solve(A.rhs().col(0));
    slv_time     += timer.stop();
    v_sol.extract(Psi);
    for(index_t Mp=0; Mp<mp.dim(); ++Mp){
    gsMultiPatch<> PsiPsitp_temp;
    timer.restart();
    vsolVector    = solver.solve(A.rhs().col(Mp));
    slv_time     += timer.stop();
    v_sol.extract(PsiPsitp_temp);
    Psi.embed(Mp+1);
    Psi.patch(0).coefs().col(Mp) = PsiPsitp_temp.patch(0).coefs().col(0);
    }
    // // ... correct boundary
    ProjectionNormalCPoints(Psi);
    Psi.addAutoBoundaries();
    Psi.computeTopology();
    auto  sv0       = solVector; // initial guess set as last solution
    solution u_lsol = A.getSolution(u, sv0);
    //! [Solver loop]
    gsInfo<< A.numDofs() <<std::flush;
    for(int ip{0}; ip<=this->m_maxIter; ++ip)
    {
        // ...  0  dirichlet for boundaries
        //u.setup(bc_mae, dirichlet::l2Projection, 0);
        // Initialize the system
        A.initSystem();
        // .. update Coeffeicient of conductivity
        auto  ExprMAE     = pow( abs(pow(div(PP).val(),IGdim) - gammaMAE*jac(PP).det())+ gammaMAE*CoeffDensity/(int_uh_0*abs(rho.val()) + int_uh_1), 1./IGdim);
        auto IntegDensity = ev.integral(ExprMAE);
        CoeffConductivity = Neumann_Int/IntegDensity;
        // MAE system
        A.assemble(u * CoeffConductivity * (-1.) * ExprMAE);//rhs vector

        // Compute the Neumann terms defined on physical space
        auto g_N = A.getBdrFunction(G);
        A.assembleBdr(bc_mae.get("Neumann"), u * g_N.tr() * nv(G) );

        gsInfo<< " ." <<std::flush;// Assemblying done
        timer.restart();
        // solVector = Poisson.solve(A.rhs());
        // solVector = Ssolver.solve(A.rhs());
        solVector = Ssolver.solveWithGuess(A.rhs(), sv0);
        slv_time += timer.stop();
        gsInfo<< "." <<std::flush; // Linear solving done

        // omp_set_dynamic(0);     // Explicitly disable dynamic teams
        // omp_set_num_threads(1); // Use these threads for later parallel regions
        auto h1Res= math::sqrt(ev.integral( ( grad(u_sol)-grad(u_lsol)).sqNorm()  ));
        sv0       = solVector;

        if ( h1Res < 1e-8 || ip == this->m_maxIter ){
            // ! end Picard loop
            gsInfo<< "\n Niter in Picard : " << ip
                    << ".. H1 residual : "<<std::scientific<<h1Res
                    << ".. min JAcobian : "<<ev.min(jac(PP).det())<<"..";
            break;
            } //
        A.initSystem(IGdim);
        // Obtain control points for the gradient of Psi
        A.assemble(u * grad(u_sol));
        for(index_t Mp=0; Mp<mp.dim(); ++Mp){
        gsMultiPatch<> PsiPsitp_temp;
        timer.restart();
        // vsolVector = Poisson.L2ProjectScalar(A.rhs().col(Mp));
        // vsolVector = solver.solve(A.rhs().col(Mp));
        vsolVector    = solver.solveWithGuess(A.rhs().col(Mp), Psi.patch(0).coefs().col(Mp));
        slv_time     += timer.stop();
        v_sol.extract(PsiPsitp_temp);
        Psi.patch(0).coefs().col(Mp) = PsiPsitp_temp.patch(0).coefs().col(0);
        }
        // // ... correct boundary
        ProjectionNormalCPoints(Psi);
        Psi.addAutoBoundaries();
        Psi.computeTopology();
    }//for loop
    // ...
    timer.stop();
    gsInfo<<" CPU-time : "<< slv_time   <<"<>\n";
    return Psi;

};

// computes the projection L^2 of a composition and return a MultiPatch object
gsMultiPatch<> gsAdaptiveMultiPatchBuilder::buildCompMultiPatch(gsMultiPatch<> Psitp, index_t elevDegree, double quadValue) const 
{
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;
    gsSparseSolver<>::CGDiagonal solver;

    gsInfo<<"<> computes composition \n";
    gsMultiPatch<> Psi;
    
    double slv_time(0);
    gsStopwatch timer;
    timer.restart();
    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    // It could be beneficial for the composition of the two mappings
    A.options().setReal("quA", quadValue);
    // A.options().setInt("quRule", 2);
    A.options().setSwitch("SameElement",false); // Very important for the composition of the two mappings
    // Elements used for numerical integration
    gsMultiBasis<> basis_comp = this->n_basis;
    basis_comp.degreeElevate(elevDegree);
    A.setIntegrationElements(basis_comp);

    //... 
    space v        = A.getSpace(basis_comp);
    gsMatrix<> vsolVector;
    solution v_sol = A.getSolution(v, vsolVector);

    //::::::::::::::::::::    Compute the composition of geometry maps      :::::::::::::::::::::::::
    geometryMap PP = A.getMap(Psitp);
    //...
    A.initSystem();
    A.assemble(v*v.tr());//Matrix in one dimension
    solver.compute( A.matrix() );

    auto comp = A.getCoeff(this->m_mapping, PP);

    A.initSystem(this->m_mapping.geoDim());
    //Obtain control points for the gradient of mpLeft.comp(Psi)
    A.assemble(v * comp.tr() );// blocked by this one
    vsolVector = solver.solve(A.rhs().col(0));
    v_sol.extract(Psi);

    for(index_t i=1; i<this->m_mapping.geoDim(); ++i)
    {
        gsMultiPatch<> PsiVec;
        vsolVector = solver.solve(A.rhs().col(i));
        v_sol.extract(PsiVec);
        Psi.embed(i+1);
        Psi.patch(0).coefs().col(i) = PsiVec.patch(0).coefs();
    }
    Psi.addAutoBoundaries();
    Psi.computeTopology();
    //#-++++++++++++++++++++++++ End of sharing part of any geometry------------------------------
    slv_time += timer.stop();
    timer.stop();
    gsInfo<<" CPU-time : "<< slv_time   <<"<>\n";
    return Psi;
};


// computes the projection L^2 of a composition and return a MultiPatch object
gsMultiPatch<> gsAdaptiveMultiPatchBuilder::buildCompBasisMultiPatch(const gsMultiBasis<> dbasis, const gsMultiPatch<> Psitp) const 
{
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;
    gsSparseSolver<>::CGDiagonal solver;

    gsInfo<<"<> computes composition \n";
    gsMultiPatch<> Psi;
    
    double slv_time(0);
    gsStopwatch timer;
    timer.restart();
    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    // A.options().setInt("quRule", 2);
    A.options().setSwitch("SameElement",false); // Very important for the composition of the two mappings

    A.setIntegrationElements(dbasis);
    //... 
    space v        = A.getSpace(dbasis);
    gsMatrix<> vsolVector;
    solution v_sol = A.getSolution(v, vsolVector);

    //::::::::::::::::::::    Compute the composition of geometry maps      :::::::::::::::::::::::::
    geometryMap PP = A.getMap(Psitp);
    //...
    A.initSystem();
    A.assemble(v*v.tr());//Matrix in one dimension
    solver.compute( A.matrix() );

    auto comp = A.getCoeff(this->m_mapping, PP);

    A.initSystem(this->m_mapping.geoDim());
    //Obtain control points for the gradient of mpLeft.comp(Psi)
    A.assemble(v * comp.tr() );// blocked by this one
    vsolVector = solver.solve(A.rhs().col(0));
    v_sol.extract(Psi);

    for(index_t i=1; i<this->m_mapping.geoDim(); ++i)
    {
        gsMultiPatch<> PsiVec;
        vsolVector = solver.solve(A.rhs().col(i));
        v_sol.extract(PsiVec);
        Psi.embed(i+1);
        Psi.patch(0).coefs().col(i) = PsiVec.patch(0).coefs();
    }
    Psi.addAutoBoundaries();
    Psi.computeTopology();
    //#-++++++++++++++++++++++++ End of sharing part of any geometry------------------------------
    slv_time += timer.stop();
    timer.stop();
    gsInfo<<" CPU-time : "<< slv_time   <<"<>\n";
    return Psi;
};

//.......................................................
//........... useful functions for time moving meshes
//.......................................................


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

    // dbasis0.assign(degree + 1, 0.0);
    // dbasis1.assign(degree + 1, 0.0);
    //dbasis2.assign(degree + 1, 0.0); // second derivative not used

    index_t minm = std::min(degree,2);// if degree == 1, then all second derivatives are 0

    gsVector<double> left(degree), right(degree);
    gsMatrix<double> ndu(degree + 1, degree + 1);
    gsMatrix<double> a(2, degree + 1);
    gsMatrix<double> ders(3, degree + 1);
    //std::vector<std::vector<double>> ndu(degree + 1, std::vector<double>(degree + 1));
    //std::vector<std::vector<double>> a(2, std::vector<double>(degree + 1));
    //std::vector<std::vector<double>> ders(3, std::vector<double>(degree + 1, 0.0));
    span     = gsAdaptiveMultiPatchBuilder::find_span(knots, degree, x);

    ndu(0,0)        = 1.0;
    // computes basis functions
    #pragma omp parallel for
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
    #pragma omp parallel for
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
        //dbasis2(j) = ders(2,j) * degree * (degree - 1);
    }
};

// compute the right-hand side vector for the adaptive multi-patch assembly : one patch
// This function computes the right-hand side vector for a given knot span and degree
// vector_un : vector of control points for the solution
// vector_mp : vector of control points for the MAE adaptive mapping : unit square
// vector_cp : vector of control points for the adaptive multi-patch composition mapping
void gsAdaptiveMultiPatchBuilder::assemble_rhsvector_ad(const index_t& p1, const index_t& p2,
                           const gsKnotVector<double>& knots_1, const gsKnotVector<double>& knots_2,
                           const gsMatrix<double>& vector_mp, const gsMatrix<double>& vector_cp,
                           const gsMatrix<double>& vector_un, gsMatrix<double>& rhs) const {
    
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
    index_t nb12 = nb1 * nb2;
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

                            ad_x       += vector_mp(gi) * bi_0;
                            ad_y       += vector_mp(nb12+gi) * bi_0;

                            ad_xx      += vector_cp(gi) * bi_x;
                            ad_xy      += vector_cp(gi) * bj_y;
                            ad_yx      += vector_cp(nb12+gi) * bi_x;
                            ad_yy      += vector_cp(nb12+gi) * bj_y;

                            val_un     += vector_un(gi) * bi_0;
                        }
                    }
                    // basis functions in the image of quadrature points by Adaptive mapping
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

// compute the basis functions for a given knot span and degree and x value
void gsAdaptiveMultiPatchBuilder::nurbsbasis_functions(const gsKnotVector<double>& knots,const gsKnotVector<double>& omega, const index_t& degree, const double& x, index_t& span,
                     gsVector<double>& dbasis0,
                     gsVector<double>& dbasis1)  const 
{

    // dbasis0.assign(degree + 1, 0.0);
    // dbasis1.assign(degree + 1, 0.0);
    //dbasis2.assign(degree + 1, 0.0); // second derivative not used

    index_t minm = std::min(degree,2);// if degree == 1, then all second derivatives are 0

    gsVector<double> left(degree), right(degree);
    gsMatrix<double> ndu(degree + 1, degree + 1);
    gsMatrix<double> a(2, degree + 1);
    gsMatrix<double> ders(3, degree + 1);
    //std::vector<std::vector<double>> ndu(degree + 1, std::vector<double>(degree + 1));
    //std::vector<std::vector<double>> a(2, std::vector<double>(degree + 1));
    //std::vector<std::vector<double>> ders(3, std::vector<double>(degree + 1, 0.0));
    span     = gsAdaptiveMultiPatchBuilder::find_span(knots, degree, x);

    ndu(0,0)        = 1.0;
    // computes basis functions
    #pragma omp parallel for
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
    #pragma omp parallel for
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
    double sum_0 = 0.;
    double sum_x = 0.;
    for (index_t j = 0; j <= degree; ++j) {
        dbasis0(j) = ders(0,j);//*omega(span-degree+j);
        sum_0     += dbasis0(j);
        dbasis1(j) = ders(1,j) * degree;//*omega(span-degree+j);
        sum_x     += dbasis1(j);
        //dbasis2(j) = ders(2,j) * degree * (degree - 1);
    }
    // for (index_t j = 0; j <= degree; ++j) {
    //     dbasis0(j) = dbasis0(j)/sum_0;//R_j
    //     dbasis1(j) = (dbasis1(j)-dbasis0(j) * sum_x)/sum_0;//R'_j
    // }
};

// compute the right-hand side vector for the adaptive multi-patch assembly : one patch
// This function computes the right-hand side vector for a given knot span and degree
// vector_un : vector of control points for the solution
// vector_mp : vector of control points for the MAE adaptive mapping : unit square
// vector_cp : vector of control points for the adaptive multi-patch composition mapping
void gsAdaptiveMultiPatchBuilder::assemble_nurbsrhsvector_ad(const index_t& p1, const index_t& p2,
                           const gsKnotVector<double>& knots_1, const gsKnotVector<double>& knots_2,
                           const gsKnotVector<double>& omega_1, const gsKnotVector<double>& omega_2,
                           const gsMatrix<double>& vector_mp, const gsMatrix<double>& vector_cp,
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
    index_t nb12 = nb1 * nb2;

    // Compute the number of elements in each direction
    index_t ne1 = 10*(nb1 - p1);
    index_t ne2 = 10*(nb2 - p2);
    double  h1  = 0.1*knots_1[p1+1];
    double  h2  = 0.1*knots_2[p2+1];

    #pragma omp parallel for
    for (index_t ie1 = 0; ie1 < ne1; ++ie1) {

        // coefs of the element for quadrature points
        double a1   = ie1*h1;//knots_1[ie1 + p1];
        double b1   = ie1*h1+h1;//knots_1[ie1 + p1 + 1];
        double c0_1 = 0.5 * (a1 + b1);
        double c1_1 = 0.5 * (b1 - a1);

        for (index_t ie2 = 0; ie2 < ne2; ++ie2) {

            double a2   = ie2*h2;//knots_2[ie2 + p2];
            double b2   = ie2*h2+h2;//knots_2[ie2 + p2 + 1];
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
                nurbsbasis_functions(knots_1, omega_1, p1, x1, span1, xbasis_0, xbasis_1);

                for (index_t g2 = 0; g2 < m; ++g2) {
                    // map the quadrature points to the element
                    double x2 = c1_2 * u(g2) + c0_2;
                    double w2 = c1_2 * w(g2);
                    // Compute the basis functions in the second direction
                    index_t span2;
                    gsVector<double> ybasis_0(p2 + 1);
                    gsVector<double> ybasis_1(p2 + 1);
                    nurbsbasis_functions(knots_2, omega_2, p2, x2, span2, ybasis_0, ybasis_1);

                    // Assembles adaptive mesh
                    double ad_x   = 0.0;
                    double ad_y   = 0.0;
                    for (index_t i = 0; i <= p1; ++i) {
                        for (index_t j = 0; j <= p2; ++j) {
                            index_t gi  = (span1 - p1+i) + nb1*(span2 - p2+j);

                            double bi_0 = xbasis_0(i) * ybasis_0(j);

                            ad_x       += vector_cp(gi) * bi_0;
                            ad_y       += vector_cp(nb12+gi) * bi_0;
                        }
                    }
                    // basis functions in the image of quadrature points by Adaptive mapping
                    index_t ad_span1;
                    gsVector<double> ad_xbasis_0(p1 + 1);
                    gsVector<double> ad_xbasis_1(p1 + 1);
                    nurbsbasis_functions(knots_1, omega_2, p1, ad_x, ad_span1, ad_xbasis_0, ad_xbasis_1);
                    index_t ad_span2;
                    gsVector<double> ad_ybasis_0(p2 + 1);
                    gsVector<double> ad_ybasis_1(p2 + 1);
                    nurbsbasis_functions(knots_2, omega_2, p2, ad_y, ad_span2, ad_ybasis_0, ad_ybasis_1);
                    // Assembles adaptive mesh
                    ad_x   = 0.0;
                    ad_y   = 0.0;
                    for (index_t i = 0; i <= p1; ++i) {
                        for (index_t j = 0; j <= p2; ++j) {
                            index_t gi  = (ad_span1 - p1+i )+ nb1*(ad_span2 - p2+j);

                            double bi_0 = ad_xbasis_0(i) * ad_ybasis_0(j);

                            ad_x       += vector_mp(gi) * bi_0;
                            ad_y       += vector_mp(nb12+gi) * bi_0;
                        }
                    }

                    double weight = w1 * w2;
                    for (index_t i = 0; i <= p1; ++i) {
                        for (index_t j = 0; j <= p2; ++j) {
                            index_t gi   = (span1 - p1+i )+ nb1*(span2 - p2+j);
                            double bi_0  = xbasis_0(i) * ybasis_0(j);

                            rhs(gi)      += ad_x * bi_0 * weight;
                            rhs(nb12+gi) += ad_y * bi_0 * weight;
                        }
                    }
                }
            }
        }
    }
};