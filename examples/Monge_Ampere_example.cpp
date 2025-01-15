/** @file Monge_Ampere_example.cpp

    @brief Tutorial on how to use expression assembler to solve a non-linear Monge-Ampere equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris & M. BAHARI
*/

//! [Include namespace]
#include <gismo.h>

using namespace gismo;
//! [Include namespace]

template<typename T>
std::vector< gsSparseMatrix<T> > assembleTensorMass(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt
    )
{
    const index_t d = basis.dim();
    std::vector< gsSparseMatrix<T> > result(d);
    for ( index_t i=0; i!=d; ++i )
    {
        result[i] = gsPatchPreconditionersCreator<>::massMatrix(basis.component(d-1-i));
        //eliminateDirichlet1D(boundaryConditionsForDirection(bc,d-1-i), opt, result[i]);
    }
    return result;
}

template<typename T>
std::vector< gsSparseMatrix<T> > assembleTensorStiffness(
    const gsBasis<T>& basis,
    const gsBoundaryConditions<T>& bc,
    const gsOptionList& opt
    )
{
    const index_t d = basis.dim();
    std::vector< gsSparseMatrix<T> > result(d);
    for ( index_t i=0; i!=d; ++i )
    {
        result[i] = gsPatchPreconditionersCreator<>::stiffnessMatrix(basis.component(d-1-i));
        //eliminateDirichlet1D(boundaryConditionsForDirection(bc,d-1-i), opt, result[i]);
    }
    return result;
}

template<typename T>
class Poisson_FastDiag {
public:
    // This class implements the fast diagonalization algorithm 
    // described in the paper: https://doi.org/10.1016/j.cma.2023.116570
    Poisson_FastDiag(const gsBasis<T>& basis,
            const gsBoundaryConditions<T>& bc,
            const gsOptionList& opt,
            T tau = 0.0)
        : _tau(tau)
    {
        const index_t rdim = basis.dim();

        // Assemble univariate local stiff matrix
        std::vector< gsSparseMatrix<T> > Ks = assembleTensorStiffness(basis, bc, opt);
        // Assemble univariate local mass matrix
        std::vector< gsSparseMatrix<T> > Ms  = assembleTensorMass(basis, bc, opt);

        typedef typename gsMatrix<T>::GenSelfAdjEigenSolver EVSolver;
        EVSolver es;

        for (size_t i = 0; i < Ms.size(); ++i) {
            /*
            We first consider the generalized eigendecompositions problems
            𝐾1 𝑈1 = 𝑀1𝑈1d1, 𝐾2𝑈2 = 𝑀2𝑈2d2, 𝐾3𝑈3 = 𝑀3𝑈3d3,
            where d1, d2 and d3 are diagonal matrices such that 𝑈1, 𝑈2 and 𝑈3 fulfills
            𝑈^𝑇𝑀_1𝑈_1 = 𝐼1, 𝑈^𝑇𝑀_2𝑈_2 = 𝐼_2, 𝑈^𝑇𝑀_3𝑈_3 = 𝐼_3
            */
            //gsEigen::GeneralizedSelfAdjointEigenSolver<gsSparseMatrix<T>> es(Ks[i], Ms[i]);
            es.compute(Ks[i], Ms[i], gsEigen::ComputeEigenvectors);

            ds.push_back(es.eigenvalues());
            Us.push_back(es.eigenvectors());
            t_Us.push_back(es.eigenvectors().transpose());
        }
        // TODO :  avoid kron !
        if (rdim == 2) {
            forward  = t_Us[0].kron(t_Us[1]).eval();
            backward = Us[0].kron(Us[1]).eval();
        } else {
            forward  = t_Us[0].kron(t_Us[1].kron(t_Us[2])).eval();
            backward = Us[0].kron(Us[1].kron(Us[2])).eval();
        }
        
        M_proj = backward * forward;        
        this->_rdim = rdim;
    }
    mutable gsMatrix<T> s_tilde;
    mutable gsMatrix<T> r_tilde;

    gsMatrix<T> solve(const gsMatrix<T>& b) const {
        if (_rdim == 2) {
            index_t n1 = ds[0].rows();
            index_t n2 = ds[1].rows();

            s_tilde = b.reshape(n1,n2);
            s_tilde = t_Us[1] * s_tilde *Us[0];
            #pragma omp parallel for
            for (index_t i1 = 0; i1 < n1; ++i1) {
                for (index_t i2 = 0; i2 < n2; ++i2) {
                    s_tilde(i1, i2) = s_tilde(i1, i2) / (ds[0](i1,0) + ds[1](i2,0) + _tau);
                }
            }
            s_tilde = Us[1] * s_tilde * t_Us[0];
            s_tilde = s_tilde.reshape(n1*n2, 1);
            return s_tilde;
        } else {
            index_t n1 = ds[0].rows();
            index_t n2 = ds[1].rows();
            index_t n3 = ds[3].rows();

            s_tilde = b.reshape(n1,n2*n3);
            s_tilde = t_Us[2] * s_tilde * Us[1]* Us[0];
            #pragma omp parallel for
            for (index_t i1 = 0; i1 < n1; ++i1) {
                for (index_t i2 = 0; i2 < n2; ++i2) {
                    for (index_t i3 = 0; i3 < n3; ++i3) {
                        index_t k = i3 + i2 * n3;
                        s_tilde(i1,k) = s_tilde(i1, k) / (ds[0](i1,0) + ds[1](i2,0) + ds[2](i3,0) + _tau);
                    }
                }
            }            
            s_tilde = Us[2] * s_tilde * t_Us[1]* t_Us[0];
            s_tilde = s_tilde.reshape(n1*n2*n3, 1);
            return s_tilde;
        }

    }
    // Computes the L2-projection of a function using M_proj * ∫(funct * v),
    // where M_proj is the mass matrix inverse, funct is the input function, and v is the test function.
    gsMatrix<T> L2ProjectScalar(const gsMatrix<T>& b) const {
        if(_rdim == 2){
            index_t n1 = ds[0].rows();
            index_t n2 = ds[1].rows();

            s_tilde = b.reshape(n1,n2);
            s_tilde = t_Us[1] * s_tilde *Us[0];

            s_tilde = Us[1] * s_tilde * t_Us[0];
            s_tilde = s_tilde.reshape(n1*n2, 1);
            return s_tilde;
        } else{
            index_t n1 = ds[0].rows();
            index_t n2 = ds[1].rows();
            index_t n3 = ds[3].rows();

            s_tilde = b.reshape(n1,n2*n3);
            s_tilde = t_Us[2] * s_tilde * Us[1]* Us[0];

            s_tilde = Us[2] * s_tilde * t_Us[1]* t_Us[0];
            s_tilde = s_tilde.reshape(n1*n2*n3, 1);
            return s_tilde;
        }
    }
    // Computes the L2-projection of a function using M_proj * ∫(funct * v),
    // where M_proj is the mass matrix inverse, funct is the input function, and v is the test function.
    gsMatrix<T> L2ProjectVec(const gsMatrix<T>& b, bool other = false) const {
        // other = true : the surface is in 3D, set _rdim to 3.
        s_tilde = b;
        index_t nsize = (int)(b.size()/_rdim);
        if(other){
            nsize = (int)(b.size()/3);
            s_tilde.reshape(nsize,3);
            s_tilde.col(0) = M_proj *s_tilde.col(0);
            s_tilde.col(1) = M_proj *s_tilde.col(1); 
            s_tilde.col(2) = M_proj *s_tilde.col(2); 
            s_tilde.reshape(nsize*2,1);
            return s_tilde;
        }
        if (_rdim == 2) {
            // TODO

            s_tilde.reshape(nsize,2);
            s_tilde.col(0) = M_proj *s_tilde.col(0);
            s_tilde.col(1) = M_proj *s_tilde.col(1); 
            s_tilde.reshape(nsize*2,1);


        } else {
            // TODO
            s_tilde.reshape(nsize,3);
            s_tilde.col(0) = M_proj *s_tilde.col(0);
            s_tilde.col(1) = M_proj *s_tilde.col(1); 
            s_tilde.col(2) = M_proj *s_tilde.col(2); 
            s_tilde.reshape(nsize*3,1);
        }
        return s_tilde;
    }
private:

    std::vector<gsMatrix<T>> ds;
    std::vector<gsMatrix<T>> Us, t_Us;
    gsMatrix<T> M_proj, forward, backward;
    int _rdim;
    T _tau;
};

void ProjectionNormalCPoints(gsMultiPatch<>& Psi, int boxMaxNumber = 1){
    // Projection normal of control points (exact geometry)
    for (int boxNumber = 0; boxNumber < boxMaxNumber; ++boxNumber)
    {
        // test if the boundary interface is not an inner interface between patches
        auto lVal = int(1.1*Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]);
        auto hVal = int(1.1*Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(2).at(0) ).array()[0]);
        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(1).size(); ++i_x) // x=0 control points be like (0,:) in this case
        {
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[0] = lVal;
        }

        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(2).size(); ++i_x)// x=1 control points be like (1,:) in this case
        {
        Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(2).at(i_x) ).array()[0] = hVal;
        }

        lVal = int(1.1*Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(3).at(0) ).array()[1]);
        hVal = int(1.1*Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(4).at(0) ).array()[1]);
        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(3).size(); ++i_x) // y=0 control points be like (:,0) in this case
        {
        Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(3).at(i_x) ).array()[1] = lVal;
        }
        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(4).size(); ++i_x)// y=1 control points be like (:,1) in this case
        {
        Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(4).at(i_x) ).array()[1] = hVal;
        }
    }
};


int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot           = false;
    index_t numRefine   = 3;
    index_t numLRefine  = 3;
    index_t numElevate  = 1;
    index_t maxIter     = 30;
    double eps          = 1e-5; // pinalization coefficient
    double tolPicard    = 1e-8;
    double IntensityMAE = 10.;
    bool plotMAeRes     = false;
    bool export_b64     = false;
    // Specify the file path
    //std::string fn("pde/quart_annulus.xml");
    //std::string fn("pde/infinit_plate.xml");
    std::string fn("pde/circle.xml");

    gsCmdLine cmd("Tutorial on solving a non-linear Monge-Ampere problem.");
    cmd.addInt("i", "iter", "Maximum number of iterations for the iterative Picard", maxIter);
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "u", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addInt( "l", "numLRefine", "Number of local h-refinement loops",  numLRefine );
    cmd.addString( "d", "file", "Input XML file data", fn );
    cmd.addInt("quRule",
                 "Quadrature rule [1:GaussLegendre,2:GaussLobatto,3:PatchRule]",
                 1);
    cmd.addSwitch("plotMAeRes", "PLot only result of solving MA equation", plotMAeRes);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addReal( "f", "IntensityMAE", "Intensity of density function",  IntensityMAE);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Load the file
    gsFileData<> fd(fn);
    gsInfo << "Loaded file " << fd.lastPath() << "\n";
    // Create a gsMultipatch and add the loaded geometry
    gsMultiPatch<> mpLeft; // = gsNurbsCreator<>::BSplineSquareGrid(1,1,1, 0.0, 0.0);
    fd.getId(1,mpLeft);
    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    mpLeft.degreeElevate(numElevate);
    mpLeft.computeTopology();

    // .... one single patch
    //gsInfo << "INFO IN PARAMETRIC DOMAIN "<< mpLeft.dim() << mpLeft.parDim() <<"\n";
    gsMultiPatch<> mp = gsNurbsCreator<>::BSplineSquareGrid(1,1,1, 0.0, 0.0);
    //Get all interfaces and boundaries:
    mp.degreeElevate(numElevate);
    mp.computeTopology();
    //mp.addAutoBoundaries();

    //..... Test 2
    // Manufactured identity mapping
    gsFunctionExpr<> sN("x","y",2);
    // Right-hand side function : Analytical density function (det(H(u))=f= sigma/rho)
    // Load the file
    gsFunctionExpr<> f;
    fd.getId(2003, f);
    gsInfo<<"Density function "<< f << "\n";

    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);
    // For simplicity, set Neumann boundary conditions
   for ( gsMultiPatch<>::const_biterator
            bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
   {
       bc.addCondition( *bit, condition_type::neumann, &sN );
   }
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";

    //! [Refinement]
    gsMultiBasis<double> dbasis(mpLeft, true);//true: poly-splines (not NURBS)
    
    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
#ifdef _OPENMP
    gsInfo<< "Available threads: "<< omp_get_max_threads() <<"\n";
#endif
    //! [Refinement]
    
    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    //A.setOptions(Aopt);
    gsInfo<<"Active options:\n"<< A.options() <<"\n";

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);
    // Set the geometry map
    geometryMap G = A.getMap(mp);

    // Set the Target geometry map
    geometryMap GLeft = A.getMap(mpLeft);

    // Set pow for BFO method dim in parameteric domain
    auto IGdim     = G.domainDim();

    // Set dimension of target geometry
    auto ITdim     = mpLeft.geoDim();

    // Set factor for BFO method
    auto gammaMAE = factorial(G.domainDim());

    // Set the discretization space
    space u = A.getSpace(dbasis);
    
    // Set the source term with respect to target geometry
    auto ff = A.getCoeff(f, GLeft);

    //gsFunctionExpr<> sI("0.5*(x**2+y**2)+x*y",2);
    auto u_I = ev.getVariable(sN, G);

    // Solution vector and solution variable
    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);
    //! [Problem setup]

    //! [Solver loop]
    gsSparseSolver<>::CGDiagonal solver;
    
    gsInfo<< "(dot1=assembled, dot2=solved)\n"
        "\nDoFs: ";
    double setup_time(0), ma_time(0), slv_time(0);    
    gsStopwatch timer;
    timer.restart();
    //::::::::::::::::::::      mesh adaptation solver         :::::::::::::::::::::::::
    for (int r=0; r<=numRefine; ++r)
    {
        dbasis.uniformRefine();
        mp.uniformRefine();
        mpLeft.uniformRefine();
    }
    
    //auto Poisson  = gsPatchPreconditionersCreator<>::fastDiagonalizationOp(dbasis.basis(0),bc,A.options(), 1.,eps,0.);
    Poisson_FastDiag<double> Poisson(dbasis.basis(0), bc, A.options(), eps);

    u.setup(bc, dirichlet::l2Projection, 0);
    // Compute the system matrix and right-hand side

    // Initialize the system :  identity mapping as initial guess
    A.initSystem();
    // Initialize the system : start Computing the conductivity coeffeicient ...
    // Compute the Neumann terms defined on physical space
    gsInfo << "evaluate integral " << ev.integral(ff.val()) << "\n";

    auto g_N = A.getBdrFunction(G);
    auto Neumann_Int{ev.integralBdrBc(bc.get("Neumann"), g_N.tr() * nv(G) )};
    //... nromalisation of density function
    auto CoeffDensity{ev.integral((1.+IntensityMAE*ff.val()))};
    auto CoeffConductivity{Neumann_Int/ev.integral(pow(pow(IGdim,IGdim)+gammaMAE * CoeffDensity/(1.+IntensityMAE*ff.val()), 1./IGdim) )};

    setup_time += timer.stop();

    gsInfo<< A.numDofs() <<std::flush;

    timer.restart();
    A.assemble(
    igrad(u, G) * igrad(u, G).tr()  + eps * u *u.tr() //matrix
    ,
    u*  CoeffConductivity * (-1.)*pow(pow(IGdim,IGdim)+gammaMAE* CoeffDensity/(1.+IntensityMAE*ff.val()), 1./IGdim)  //rhs vector
    );
    
    // Compute the Neumann terms defined on physical space
    //auto g_N = A.getBdrFunction(G);
    A.assembleBdr(bc.get("Neumann"), u * g_N.tr() * nv(G) );
    A.assembleIfc(mp.interfaces(), u.left() * (u_I.tr() * nv(G.left())));
    A.assembleIfc(mp.interfaces(), u.right() * (u_I.tr() * nv(G.right())));

    ma_time += timer.stop();

    gsInfo<< "." <<std::flush;// Assemblying done
    solVector = A.rhs();
    timer.restart();
    solVector = Poisson.solve(A.rhs());

    // solver.compute( A.matrix() );
    // solVector = solver.solve(A.rhs());
    slv_time += timer.stop();

    gsInfo<< "." <<std::flush; // Linear solving done

    // Picard loop
    index_t NiterPicard{0};
    gsMatrix<> sv0; //
    solution u_lsol = A.getSolution(u, sv0);
    for(int ip{0}; ip<=maxIter; ++ip)
    {
        timer.restart();
        gsMultiPatch<> UU;
        u_sol.extract(UU);
        gsWrite(UU, "U_solution");
        auto u_s = A.getCoeff(UU);

        space v = A.getSpace(dbasis);
        gsMatrix<> vsolVector;
        solution v_sol = A.getSolution(v, vsolVector);
        A.initSystem(2);

        // Obtain control points for the gradient of Psi
        A.assemble( v * v.tr() , v * igrad(u_s,G) );
        //gsInfo <<"rhs vec = " << A.rhs().size() << "\n";
        vsolVector = Poisson.L2ProjectVec(A.rhs());
        
        gsMultiPatch<> Psi, PsiLoc;
        v_sol.extract(Psi);
        v_sol.extract(PsiLoc);

        // ... correct boundary
        ProjectionNormalCPoints(Psi);
        Psi.addAutoBoundaries();
        Psi.computeTopology();
        geometryMap PP    = A.getMap(Psi);
        geometryMap PPLoc = A.getMap(PsiLoc);
        //::::::::::::::::::::    Compute the composition of geometry maps      :::::::::::::::::::::::::
        auto  comp = PPLoc(mpLeft);
        A.initSystem(ITdim);
        //Obtain control points for the gradient of mpLeft.comp(Psi)
        A.assemble( v * v.tr() , v * comp.tr() );// blocked by this one
        vsolVector = Poisson.L2ProjectVec(A.rhs());
        //vsolVector = solver.compute(A.matrix()).solve(A.rhs());
        v_sol.extract(PsiLoc);
        //::::::::::::::::::::      end       ::::::::::::::::::::::::: 
        geometryMap PPfLoc = A.getMap(PsiLoc);
        auto ff = A.getCoeff(f, PPfLoc);

        // auto ff = A.getCoeff(f, GLeft, PP);
        // auto ffG = A.getCoeff(f, GLeft);
        // gsInfo << "ff" << ev.integral(ff.val()) << "\n";
        // gsInfo << "ffG" << ev.integral(ffG.val()) << "\n";

        // ...  0  dirichlet for boundaries
        sv0 = solVector;
        u.setup(bc, dirichlet::l2Projection, 0);
    
        solution u_sol = A.getSolution(u, solVector);

        // Initialize the system
        A.initSystem();
        setup_time += timer.stop();

        //gsInfo<< A.numDofs() <<std::flush;

        timer.restart();
        // Compute the system matrix and right-hand side ... Monge-Ampere eqaution .....
        
        // .. update Coeffeicient of conductivity
        auto  ExprMAE = pow( pow(div(PP).val(),IGdim) + gammaMAE*(CoeffDensity/(1.+IntensityMAE*ff.val()) - jac(PP).det()), 1./IGdim);
        auto IntegDensity = ev.integral(ExprMAE);
        CoeffConductivity = Neumann_Int/IntegDensity;
        // MAE system
        //gsInfo << " end value coeff "<< CoeffConductivity << "\n";
        A.assemble(
        igrad(u, G) * igrad(u, G).tr()  +  eps * u * u.tr()//matrix
        ,
        u * CoeffConductivity * (-1.) * ExprMAE  //rhs vector
        );
        //gsInfo << "End Assemnles \n";

        // Compute the Neumann terms defined on physical space
        auto g_N = A.getBdrFunction(G);
        A.assembleBdr(bc.get("Neumann"), u * g_N.tr() * nv(G) );
        A.assembleIfc(mp.interfaces(), u.left() * (u_I.tr() * nv(G.left())));
        A.assembleIfc(mp.interfaces(), u.right() * (u_I.tr() * nv(G.right())));
        ma_time += timer.stop();

        // gsDebugVar(A.matrix().toDense());
        // gsDebugVar(A.rhs().transpose()   );
        

        gsInfo<< " ." <<std::flush;// Assemblying done

        timer.restart();
        solVector = Poisson.solve(A.rhs());
        slv_time += timer.stop();

        gsInfo<< "." <<std::flush; // Linear solving done

        // omp_set_dynamic(0);     // Explicitly disable dynamic teams
        // omp_set_num_threads(1); // Use these threads for later parallel regions

        ++NiterPicard;
        auto l2errRes = math::sqrt(ev.integral( ( igrad(u_lsol,G) - igrad(u_sol,G) ).sqNorm()  ));
        if ( l2errRes < tolPicard || ip == maxIter ){
            // ! end Picard loop
            gsInfo<< "\n Niter in Picard : " << NiterPicard << ".. L2 residual : "<<std::scientific<<l2errRes<<"\n";
            break; 
            } // 
    }//for loop
        // omp_set_dynamic(0);     // Explicitly disable dynamic teams
        // omp_set_num_threads(1); // Use these threads for later parallel regions
 
    timer.stop();
    gsInfo<<"\n\nTotal time: "<< setup_time+ma_time+slv_time<<"\n";
    gsInfo<<"     Setup: "<< setup_time <<"\n";
    gsInfo<<"  Assembly: "<< ma_time    <<"\n";
    gsInfo<<"   Solving: "<< slv_time   <<"\n";

    if(plotMAeRes){
        gsInfo<<"Plotting in Paraview...\n";
        gsParaviewCollection collection("ParaviewOutput/solution", &ev);
        collection.options().setSwitch("plotElements", true);
        collection.options().setSwitch("base64", export_b64);
        collection.options().setInt("plotElements.resolution", 16);
        //collection.options().setInt("numPoints", 1000);
        collection.newTimeStep(&mp);
        collection.addField(u_sol,"numerical solution");
        collection.addField(igrad(u_sol,G),"gradient_numerical solution");
        collection.addField(ff, "density function");
        collection.addField(ihess(u_sol,G).det(), "Jacobian function");
        if(maxIter == 0)
        collection.addField(CoeffConductivity * (-1.)*pow(pow(IGdim,IGdim)+gammaMAE * CoeffDensity/(1.+IntensityMAE*ff.val()), 1./IGdim) , "MAE_rhs");
        else
        collection.addField(CoeffConductivity * (-1.) * pow( pow(ilapl(u_sol,G).val(),IGdim) + gammaMAE*(CoeffDensity/(1.+IntensityMAE*ff.val()) - ihess(u_sol,G).det()), 1./IGdim) , "MAE_rhs");
        collection.saveTimeStep();
        collection.save();
        gsFileManager::open("ParaviewOutput/solution.pvd");
        }
    //! [Export visualization in ParaView]
    if (plot)
    {
        gsMultiPatch<> UU;
        u_sol.extract(UU);
        gsWrite(UU, "U_solution");
        auto u_s = A.getCoeff(UU);

        //gsMultiBasis<> gbasis(dbasis);
        //gbasis.reduceContinuity(1);
        space v = A.getSpace(dbasis);
        gsMatrix<> vsolVector;
        solution v_sol = A.getSolution(v, vsolVector);
        A.initSystem(IGdim);

        //gsVector<> pt(2); pt.setConstant(0.5);
        //ev.testEval( v, pt );
        //ev.testEval( igrad(u_sol,G), pt );

        // Obtain control points for the gradient of Psi
        A.assemble( v * v.tr() , v * igrad(u_s,G));
        vsolVector = Poisson.L2ProjectVec(A.rhs());
        gsMultiPatch<> Psi, Psitp;
        v_sol.extract(Psitp);
        //... correct the boundary
        ProjectionNormalCPoints(Psitp);

        //::::::::::::::::::::    Compute the composition of geometry maps      :::::::::::::::::::::::::
        // Psi.addAutoBoundaries();
        geometryMap PP = A.getMap(Psitp);
        auto  comp = PP(mpLeft);
        A.initSystem(ITdim);
        //Obtain control points for the gradient of mpLeft.comp(Psi)
        A.assemble( v * v.tr() , v * comp.tr() );// blocked by this one
        vsolVector = Poisson.L2ProjectVec(A.rhs());
        v_sol.extract(Psitp);
        Psitp.addAutoBoundaries();
        Psitp.computeTopology();
        gsInfo << "end of adaptive mapping computation\n" << Psitp<< "\n";

        for(size_t i =0; i<Psitp.nPatches(); ++i)
            Psi.addPatch(gsTHBSpline<2>( dynamic_cast<const gsTensorBSpline<2>&>(Psitp.patch(i)) ));
        Psi.addAutoBoundaries();
        Psi.computeTopology();

        //Psi.uniformRefine();
        gsMultiBasis<> dbasis(Psi, true);//true: poly-splines (not NURBS)

        geometryMap PPF = A.getMap(Psi);
        auto ff_TG      = A.getCoeff(f, PPF);
        // --------------- adaptive refinement ---------------
        // Specify cell-marking strategy...
        MarkingStrategy adaptRefCrit = PUCA;
        //MarkingStrategy adaptRefCrit = GARU;
        //MarkingStrategy adaptRefCrit = errorFraction;
        real_t adaptRefParam = 0.7;
        // Elements used for numerical integration
        A.setIntegrationElements(dbasis);
        gsExprEvaluator<> ev(A);

        for (int r=0; r<=numLRefine; ++r)
        {
        //! [beginRefLoop]
            gsInfo << "====== Loop " << r << " of "
                    <<numRefine<< " ====adapt Parameter ="<< adaptRefParam << " ======" << "\n";
            // --------------- error estimation/computation ---------------
            // Get the element-wise norms.
            ev.integralElWise( ( ff_TG ).sqNorm() );
            const std::vector<real_t> eltErrs  = ev.elementwise();
            //! [errorComputation]

            //! [adaptRefinementPart]
            // Mark elements for refinement, based on the computed local errors and
            // the refinement-criterion and -parameter.
            std::vector<bool> elMarked( eltErrs.size() );
            gsMarkElementsForRef( eltErrs, adaptRefCrit, adaptRefParam, elMarked);
            gsInfo <<"Marked "<< std::count(elMarked.begin(), elMarked.end(), true) <<" elements.\n";
            // Refine the marked elements with a 1-ring of cells around marked elements
            gsRefineMarkedElements( dbasis, elMarked, 1);
            gsRefineMarkedElements( Psi, elMarked, 1);
            }

        //::::::::::::::::::::      end       :::::::::::::::::::::::::   
        gsInfo<<"Plotting in Paraview...\n";
        gsParaviewCollection collection("ParaviewOutput/solution", &ev);
        collection.options().setSwitch("plotElements", true);
        collection.options().setSwitch("base64", export_b64);
        collection.options().setInt("plotElements.resolution", 16);
        //collection.options().setInt("numPoints", 10000);
        collection.newTimeStep(&Psi);
        collection.addField(ff_TG, "density function");
        collection.addField(jac(PPF).det(), "Jacobian function");
         collection.saveTimeStep();
        collection.save();

        gsFileManager::open("ParaviewOutput/solution.pvd");
        // gsWrite(Psi, "Psi_mapping");
        // gsInfo << "Result written in Psi_mapping.xml \n";
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;


}// end main