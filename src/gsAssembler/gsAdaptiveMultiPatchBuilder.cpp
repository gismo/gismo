#include <gismo.h>  // Include necessary GISMO headers
#include "gsAdaptiveMultiPatchBuilder.h"

// Constructor implementation
gsAdaptiveMultiPatchBuilder::gsAdaptiveMultiPatchBuilder(const gsMultiBasis<double> basis,
                            const gsMultiPatch<> mapping,
                            const index_t   numElevate,
                            index_t maxIter,
                            double IntensityMAE )
{
    gsInfo<<"\n <>r-refinement";
    this->m_basis        = basis;
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
    mp.degreeElevate(numElevate);
    mp.computeTopology();
    //mp.addAutoBoundaries();

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
    gsPatchPreconditionersCreator<double>::Poisson_FastDiag Poisson(basis.basis(0), bc_mae, A.options(), 1e-6);  
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
            auto lVal    = this->mp.patch(boxNumber).coef( this->mp.patch(boxNumber).basis().boundary(bndIter).at(0) ).array()[dim];
            for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(bndIter).size(); ++i_x) // x=0 control points be like (0,:) in this case
            {
                Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(bndIter).at(i_x) ).array()[dim] = lVal;
            }
            bndIter ++;
            auto hVal    = this->mp.patch(boxNumber).coef( this->mp.patch(boxNumber).basis().boundary(bndIter).at(0) ).array()[dim];
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
    A.assemble(
    u *u.tr() //matrix
    ,
    u* ff.val()  //rhs vector
    );
    densityVector        = this->Poisson.L2ProjectScalar(A.rhs());
    //...
    solution density_sol = A.getSolution(u, densityVector);
    gsMultiPatch<> density;
    density_sol.extract(density);
    gsInfo<<"<>\n";
    return  density;
}

// Build and return a density as a MultiPatch object from solution vector
gsMultiPatch<> gsAdaptiveMultiPatchBuilder::buildDensity(const std::vector<double> &elwiseERROR, const double eps, index_t circleN) const 
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
    if (this->m_basis.dim() == 2){
        index_t n1   = basis_0.basis(0).component(0).numElements();
        index_t n2   = basis_0.basis(0).component(1).numElements();

        for (index_t i1 = 0; i1 < n1; i1++){
            for (index_t j1 = 0; j1 < n2; j1++){
                auto i = i1*n2 + j1;
                double Elcontr = 0.;
                index_t s = 0;
                for (index_t i2 = std::max(0,i1-circleN); i2 <= std::min(n1,i1+circleN); i2++){
                    for (index_t j2 = std::max(0,j1-circleN); j2 <= std::min(n2,j1+circleN); j2++){
                        auto j   = i2*n2 + j2;
                        Elcontr += elwiseERROR[j];
                        s       += 1;
                    }
                }
                errorVector(i) = Elcontr/s;
            }
        }
    }
    else{
        index_t n1   = basis_0.basis(0).component(0).numElements();
        index_t n2   = basis_0.basis(0).component(1).numElements();
        index_t n3   = basis_0.basis(0).component(2).numElements();
        for (index_t i1 = 0; i1 < n1; i1++){
            for (index_t j1 = 0; j1 < n2; j1++){
               for (index_t k1 = 0; k1 < n3; k1++){
                    auto i         = i1*n2*n3 + j1*n3 + k1;
                    double Elcontr = 0.;
                    index_t s      = 0;
                    for (index_t i2 = std::max(0,i1-circleN); i2 <= std::min(n1,i1+circleN); i2++){
                        for (index_t j2 = std::max(0,j1-circleN); j2 <= std::min(n2,j1+circleN); j2++){
                            for(index_t k2 = std::max(0,k1-circleN); k2 <= std::min(n3,k1+circleN); k2++){
                                auto j   = i2*n2*n3 + j2*n3 + k2;
                                Elcontr += elwiseERROR[j];
                                s       += 1;
                            }
                        }
                    }
                    errorVector(i) = Elcontr/s;
               }
            }
        }
    }
    //... normalize the error vector
    const double Maxvalue   = errorVector.maxCoeff();
    const double Minvalue   = errorVector.minCoeff();
    // gsInfo << "Density function: min "<< errorVector.minCoeff() <<"/ max " << errorVector.maxCoeff() << "\n";    
    const double meanvalue  = eps*(Maxvalue + Minvalue);
    for (index_t i1 = 0; i1 < elwnumb; i1++){
        if (errorVector(i1) > Minvalue+meanvalue)
            errorVector(i1)  = Maxvalue;
        else errorVector(i1) = Minvalue;
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
    A.assemble(
    u *u.tr() //matrix
    ,
    u * error_sol  //rhs vector
    );
    densityVector        = this->Poisson.L2ProjectScalar(A.rhs());
    //...
    gsMultiPatch<> density;
    density_sol.extract(density);
    gsInfo<<"<>\n";
    return  density;
}

// Build and return a MultiPatch object
gsMultiPatch<> gsAdaptiveMultiPatchBuilder::buildMultiPatch(const gsMultiPatch<> &density, bool composition) const 
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
    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    // It could be beneficial for the composition of the two mappings
    //A.options().setReal("quA", 2.0);
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

    // Set dimension of target geometry
    auto ITdim     = this->m_mapping.geoDim();

    // Set factor for BFO method
    auto gammaMAE  = factorial(G.domainDim());

    // Solution vector and solution variable
    gsMatrix<> solVector;
    solution u_sol   = A.getSolution(u, solVector);

    // ---- manipulation of density function ----
    auto rho         = A.getCoeff(density);
    auto empldensity = (ev.max(abs(rho.val()))-ev.min(abs(rho.val())));
    double  int_uh_0 = 0.;
    double  int_uh_1 = 1.;
    if (empldensity < 1e-5|| this->m_IntensityMAE <= 1. )
    {
        //gsInfo << "Density function is constant in the domain rho = 1.\n";
    }
    else{
        int_uh_0  = (this->m_IntensityMAE-1.)/empldensity;
        int_uh_1  = (1.*ev.max(abs(rho.val()))-this->m_IntensityMAE*ev.min(abs(rho.val())))/empldensity;
        //gsInfo << "Density function is not constant in the domain\n";
    }
    gsInfo << "Density function: min "<< ev.min(int_uh_0*abs(rho.val()) + int_uh_1)<<"/ max " << ev.max(int_uh_0*abs(rho.val()) + int_uh_1) << "\n";
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

    A.assemble(
    u *u.tr() //matrix
    ,
    u*  CoeffConductivity * (-1.)*ExprMAE  //rhs vector
    );

    // Compute the Neumann terms defined on physical space
    A.assembleBdr(bc_mae.get("Neumann"), u * g_N.tr() * nv(G) );
    //A.assembleIfc(mp.interfaces(), u.left() * (u_I.tr() * nv(G.left())));
    //A.assembleIfc(mp.interfaces(), u.right() * (u_I.tr() * nv(G.right())));

    gsInfo<< "." <<std::flush;// Assemblying done
    solVector = A.rhs();
    solVector = this->Poisson.solve(A.rhs());
    gsInfo<< "." <<std::flush; // Linear solving done

    //! [Solver loop]
    gsInfo<< A.numDofs() <<std::flush;
    // Picard loop
    gsVector<>  h1Res(m_maxIter+1), l2err(m_maxIter+1), Iter_mae(m_maxIter+1);
    auto  sv0 = solVector; //
    solution u_lsol = A.getSolution(u, sv0);
    for(int ip{0}; ip<=m_maxIter; ++ip)
    {
        gsMultiPatch<> UU;
        u_sol.extract(UU);
        auto u_s       = A.getCoeff(UU);
        //space v        = A.getSpace(this->m_basis);
        gsMatrix<> vsolVector;
        solution v_sol = A.getSolution(u, vsolVector);

        A.initSystem(IGdim);
        // Obtain control points for the gradient of Psi
        A.assemble( u * u.tr() , u * grad(u_s) );
        vsolVector = this->Poisson.L2ProjectVec(A.rhs());

        gsMultiPatch<> Psi;
        v_sol.extract(Psi);
        // ... correct boundary
        ProjectionNormalCPoints(Psi);
        Psi.addAutoBoundaries();
        Psi.computeTopology();
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
        A.assemble(
        u * u.tr()//matrix
        ,
        u * CoeffConductivity * (-1.) * ExprMAE  //rhs vector
        );
        //gsInfo << "End Assemnles \n";

        // Compute the Neumann terms defined on physical space
        auto g_N = A.getBdrFunction(G);
        A.assembleBdr(bc_mae.get("Neumann"), u * g_N.tr() * nv(G) );
        //A.assembleIfc(mp.interfaces(), u.left() * (u_I.tr() * nv(G.left())));
        //A.assembleIfc(mp.interfaces(), u.right() * (u_I.tr() * nv(G.right())));

        gsInfo<< " ." <<std::flush;// Assemblying done
        solVector = this->Poisson.solve(A.rhs());
        gsInfo<< "." <<std::flush; // Linear solving done

        // omp_set_dynamic(0);     // Explicitly disable dynamic teams
        // omp_set_num_threads(1); // Use these threads for later parallel regions

        sv0               = sv0 - solVector;
        auto l2errRes = math::sqrt(ev.integral( ( grad(u_lsol)).sqNorm()  ));
        auto L2MAERes = math::sqrt(ev.integral( pow( CoeffDensity - (int_uh_0*abs(rho.val()) + int_uh_1)*jac(PP).det(),2)  ));
        auto Ddet     = ev.min(jac(PP).det());
        Iter_mae[ip]  = ip;
        h1Res[ip]     = l2errRes;// Compute the H1 residual
        l2err[ip]     = L2MAERes;// Compute the L2 error in MA equation
        if ( l2errRes < 1e-8 || ip == m_maxIter ){
            // ! end Picard loop
            gsInfo<< "\n Niter in Picard : " << ip
                    << ".. H1 residual : "<<std::scientific<<l2errRes
                    << ".. L2 MAE residual : "<<std::scientific<<L2MAERes
                    << ".. min JAcobian : "<<Ddet<<"..";
            break;
            } //
    }//for loop
    // omp_set_dynamic(0);     // Explicitly disable dynamic teams
    // omp_set_num_threads(1); // Use these threads for later parallel regions
    gsMultiPatch<> UU;
    u_sol.extract(UU);
    //...
    auto u_s       = A.getCoeff(UU);
    //... 
    space v        = A.getSpace(this->m_basis);
    gsMatrix<> vsolVector;
    solution v_sol = A.getSolution(v, vsolVector);
    A.initSystem(IGdim);
    // Obtain control points for the gradient of Psi
    A.assemble( v * v.tr() , v * grad(u_s) );
    // SOLVE ...
    vsolVector     = Poisson.L2ProjectVec(A.rhs());
    gsMultiPatch<> Psi;
    v_sol.extract(Psi);
    //... correct the boundary
    ProjectionNormalCPoints(Psi);
    Psi.addAutoBoundaries();
    Psi.computeTopology();
    if (composition){
        //::::::::::::::::::::    Compute the composition of geometry maps      :::::::::::::::::::::::::
        //... test if geometry is a surface in 3d
        bool other = false;
        if (IGdim < this->m_mapping.geoDim()){
            other = true;
        }
        // Psi.addAutoBoundaries();
        geometryMap PP = A.getMap(Psi);
        //PP(this->m_mapping);
        auto comp = A.getCoeff(this->m_mapping, PP);
        A.initSystem(ITdim);
        //Obtain control points for the gradient of mpLeft.comp(Psi)
        A.assemble( v * v.tr() , v * comp.tr() );// blocked by this one
        // vsolVector = solver.compute(A.matrix()).solve(A.rhs());
        vsolVector = this->Poisson.L2ProjectVec(A.rhs(), other);
        v_sol.extract(Psi);
        Psi.addAutoBoundaries();
        Psi.computeTopology();
        //#-++++++++++++++++++++++++ End of sharing part of any geometry------------------------------
        slv_time += timer.stop();
        timer.stop();
        gsInfo<<" CPU-time : "<< slv_time   <<"<>\n";
        return Psi;
    }
    else{
        slv_time += timer.stop();
        timer.stop();
        gsInfo<<" CPU-time : "<< slv_time   <<"<>\n";
        return Psi;
    }
};