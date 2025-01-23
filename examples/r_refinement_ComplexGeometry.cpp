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
#include <fstream> // For file operations

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
        }
        // dimension of paraemtric space
        this->_rdim = rdim;
    }
    mutable gsMatrix<T> s_tilde;
    mutable gsMatrix<T> r_tilde;
    mutable gsMatrix<T> t_tilde;

    // Computes the approximate solution of ∫-\nabla u *\nabla v + eps * u *v = ∫(funct * v) as input,
    // where funct is the input function, and v is the test function.
    gsMatrix<T> solve(const gsMatrix<T>& b) const {
        if (_rdim == 2) {
            index_t n1 = ds[0].rows();
            index_t n2 = ds[1].rows();

            s_tilde = b.reshape(n1,n2);
            s_tilde = Us[0].transpose()*s_tilde*Us[1];
            #pragma omp parallel for
            for (index_t i1 = 0; i1 < n1; ++i1) {
                for (index_t i2 = 0; i2 < n2; ++i2) {
                    s_tilde(i1, i2) = s_tilde(i1, i2) / (ds[0](i1,0) + ds[1](i2,0) + _tau);
                }
            }
            s_tilde = Us[0]*s_tilde * Us[1].transpose();
            s_tilde = s_tilde.reshape(n1*n2, 1);
            return s_tilde;
        } else {
            // TODO
            index_t n1 = ds[0].rows();
            index_t n2 = ds[1].rows();
            index_t n3 = ds[2].rows();

            s_tilde = b.reshape(n1,n2*n3);
            s_tilde = s_tilde.transpose() * Us[0];
            // matrix become (n2*n3, n1)
            #pragma omp parallel for
            for (index_t i1 = 0; i1 < n1; ++i1) {
                r_tilde         = s_tilde.col(i1);
                r_tilde         = r_tilde.reshape(n2, n3);
                r_tilde         = Us[1].transpose() * r_tilde * Us[2];
                //...
                for (index_t i2 = 0; i2 < n2; ++i2) {
                    for (index_t i3 = 0; i3 < n3; ++i3) {
                        r_tilde(i2, i3) = r_tilde(i2, i3) / (ds[0](i1,0) + ds[1](i2,0) + ds[2](i3,0) + _tau);
                    }
                }
                r_tilde         = Us[1] * r_tilde * Us[2].transpose();
                s_tilde.col(i1) = r_tilde.reshape(n2*n3,1);
            }
            s_tilde = Us[0] * s_tilde.transpose();

            s_tilde = s_tilde.reshape(n1*n2*n3, 1);
            return s_tilde;
        }
    }
    // Computes the L2-projection of a function using ∫(funct * v) as input,
    // where funct is the input function, and v is the test function.
    // (A\otimesB)x = vec(BXA^T)
    gsMatrix<T> L2ProjectScalar(const gsMatrix<T>& b) const {
        if(_rdim == 2){
            index_t n1 = ds[0].rows();
            index_t n2 = ds[1].rows();

            s_tilde = b.reshape(n1,n2);
            s_tilde = Us[0].transpose()*s_tilde*Us[1];
            //...
            s_tilde = Us[0]*s_tilde * Us[1].transpose();
            s_tilde = s_tilde.reshape(n1*n2, 1);
        } else{
            // TODO
            index_t n1 = ds[0].rows();
            index_t n2 = ds[1].rows();
            index_t n3 = ds[2].rows();

            s_tilde = b.reshape(n1,n2*n3);
            s_tilde = s_tilde.transpose() * Us[0];
            // matrix become (n2*n3, n1)
            #pragma omp parallel for
            for (index_t i1 = 0; i1 < n1; ++i1) {
                r_tilde         = s_tilde.col(i1);
                r_tilde         = r_tilde.reshape(n2, n3);
                r_tilde         = Us[1].transpose() * r_tilde * Us[2];
                //...
                r_tilde         = Us[1] * r_tilde * Us[2].transpose();
                s_tilde.col(i1) = r_tilde.reshape(n2*n3,1);
            }
            s_tilde = Us[0] * s_tilde.transpose();
            s_tilde = s_tilde.reshape(n1*n2*n3, 1);
            return s_tilde;
        }
    }
    // Computes the L2-projection of a function using ∫(funct * v) as input,
    // where funct is the input function, and v is the test function. 
    // ... replace (A\otimes B)x by vec(BXA^T)
    // ... replace (A\otimes B\otimes C)x by vec(CXBA^T)            
    gsMatrix<T> L2ProjectVec(const gsMatrix<T>& b, bool other = false) const {
        // other = true : the surface is in 3D, set _rdim to 3.
        r_tilde = b;
        if(other){
            // TODO :  mybe need kron after all
            // This is for composing surface mapping in 3D that has 3 component and square mapping
            index_t n1 = ds[0].rows();
            index_t n2 = ds[1].rows();
            index_t n3 = n2;

            // two components
            r_tilde = r_tilde.reshape(n1*n2*n3,3);
            //// ...step 1: reshape first component *****
            #pragma omp parallel for
            for (index_t i = 0; i<3; ++i){
                s_tilde = r_tilde.col(i);
                // step 2: first component            
                s_tilde = s_tilde.reshape(n1,n2*n3);
                s_tilde = s_tilde.transpose() * Us[0];
                // matrix become (n2*n3, n1)
                for (index_t i1 = 0; i1 < n1; ++i1) {
                    t_tilde         = s_tilde.col(i1);
                    t_tilde         = t_tilde.reshape(n2, n3);
                    t_tilde         = Us[1].transpose() * t_tilde * Us[2];
                    //...
                    t_tilde         = Us[1] * t_tilde * Us[2].transpose();
                    s_tilde.col(i1) = t_tilde.reshape(n2*n3,1);
                }
                s_tilde = Us[0] * s_tilde.transpose();
                // step 4: reshape
                r_tilde.col(i) = s_tilde.reshape(n1*n2*n3,1);
            }
            //... result
            r_tilde.reshape(3*n1*n2*n3, 1);
            return r_tilde;
        }
        if (_rdim == 2) {
            // 2D
            index_t n1 = ds[0].rows();
            index_t n2 = ds[1].rows();
            // two components
            r_tilde.reshape(n1*n2,2);
            //// ...step 1: reshape first component *****
            s_tilde = r_tilde.col(0);
            s_tilde = s_tilde.reshape(n1,n2);
            // step 2: first component
            s_tilde = Us[0].transpose()*s_tilde*Us[1];
            s_tilde = Us[0]*s_tilde * Us[1].transpose();
            // step 4: reshape
            r_tilde.col(0) = s_tilde.reshape(n1*n2,1);
            //// ...step 1: reshape second component *****
            s_tilde = r_tilde.col(1);
            s_tilde = s_tilde.reshape(n1,n2);
            // step 2: first component            
            s_tilde = Us[0].transpose()*s_tilde*Us[1];
            s_tilde = Us[0]*s_tilde * Us[1].transpose();
            // step 4: reshape
            r_tilde.col(1) = s_tilde.reshape(n1*n2,1);
            r_tilde.reshape(2*n1*n2,1);

        } else {
            // 3D
            index_t n1 = ds[0].rows();
            index_t n2 = ds[1].rows();
            index_t n3 = ds[2].rows();

            // two components
            r_tilde = r_tilde.reshape(n1*n2*n3,3);
            //// ...step 1: reshape first component *****
            #pragma omp parallel for
            for (index_t i = 0; i<3; ++i){
                s_tilde = r_tilde.col(i);
                // step 2: first component            
                s_tilde = s_tilde.reshape(n1,n2*n3);
                s_tilde = s_tilde.transpose() * Us[0];
                // matrix become (n2*n3, n1)
                for (index_t i1 = 0; i1 < n1; ++i1) {
                    t_tilde         = s_tilde.col(i1);
                    t_tilde         = t_tilde.reshape(n2, n3);
                    t_tilde         = Us[1].transpose() * t_tilde * Us[2];
                    //...
                    t_tilde         = Us[1] * t_tilde * Us[2].transpose();
                    s_tilde.col(i1) = t_tilde.reshape(n2*n3,1);
                }
                s_tilde = Us[0] * s_tilde.transpose();
                // step 4: reshape
                r_tilde.col(i) = s_tilde.reshape(n1*n2*n3,1);
            }
            //... result
            r_tilde.reshape(3*n1*n2*n3, 1);
        }
        return r_tilde;
    }
private:

    std::vector<gsMatrix<T>> ds;
    std::vector<gsMatrix<T>> Us;
    int _rdim;
    T _tau;
};

void ProjectionNormalCPoints(gsMultiPatch<>& Psi, int boxMaxNumber = 1){
    // Projection normal of control points (exact geometry)
    for (int boxNumber = 0; boxNumber < boxMaxNumber; ++boxNumber)
    {
        // test if the boundary interface is not an inner interface between patches
        auto lVal = 0.; //int(1.1*Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(0) ).array()[0]);
        auto hVal = 1.; //int(1.1*Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(2).at(0) ).array()[0]);
        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(1).size(); ++i_x) // x=0 control points be like (0,:) in this case
        {
            Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(1).at(i_x) ).array()[0] = lVal;
        }

        for (int i_x =0; i_x < Psi.patch(boxNumber).basis().boundary(2).size(); ++i_x)// x=1 control points be like (1,:) in this case
        {
        Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(2).at(i_x) ).array()[0] = hVal;
        }

        lVal = 0.; //int(1.1*Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(3).at(0) ).array()[1]);
        hVal = 1.; //int(1.1*Psi.patch(boxNumber).coef( Psi.patch(boxNumber).basis().boundary(4).at(0) ).array()[1]);
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
    index_t numElevate  = 0;
    index_t maxIter     = 50;
    double eps          = 1e-7; // pinalization coefficient
    double tolPicard    = 1e-8;
    double IntensityMAE = 6.;
    bool export_b64     = false;
    real_t adaptRefParam = 0.;     // ... adapt parameter.
    double FactRefPar    = 0.;    // ... adapt parameter : adaptRefParam += FactRefPar in each iter
    // ...PNormalCP: Correct the normal part of the mapping.
    bool PNormalCP      = true;
    // Specify the file path
    std::string fn("pde/quart_annulus.xml");

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
    cmd.addReal( "a", "adaptRefParam", "parameter for local h-refinement loops",  adaptRefParam );
    cmd.addReal( "p", "FactRefPar", "augement adaptRefParam with such quantity in local h-refinement loops",  FactRefPar );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addReal( "f", "IntensityMAE", "Intensity of density function",  IntensityMAE);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Load the file
    gsFileData<> fd(fn);
    gsInfo << "Loaded file " << fd.lastPath() << "\n";
    // Create a gsMultipatch and add the loaded geometry
    gsMultiPatch<> mpLeft;
    fd.getId(1,mpLeft);
    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    mpLeft.degreeElevate(numElevate);
    mpLeft.computeTopology();

    // .... one single patch
    gsMultiPatch<> mp = gsNurbsCreator<>::BSplineSquareGrid(1,1,1, 0.0, 0.0);
    //Get all interfaces and boundaries:
    mp.degreeElevate(numElevate);
    mp.computeTopology();
    //mp.addAutoBoundaries();

    //..... Test 1
    // convection coefficient:
    gsMatrix<> coeff_conv{1,2};
    // diffusion coefficient:
    //double coeff_diff = 1.;
    // reaction coefficient:
    //double coeff_reac = 0.;
    // Manufactured identity mapping
    gsFunctionExpr<> sN("x","y",2);
    // Right-hand side function : Analytical density function (det(H(u))=f= sigma/rho)
    // Load the file
    gsFunctionExpr<> f;
    fd.getId(2003, f);
    gsFunctionExpr<> s;
    fd.getId(2000, s);
    gsFunctionExpr<> rhs;
    fd.getId(2001, rhs);
    gsInfo<<"Density function "<< f << "\n";

    gsBoundaryConditions<> bcMAE;
    bcMAE.setGeoMap(mp);
    // For simplicity, set Neumann boundary conditions
   for ( gsMultiPatch<>::const_biterator
            bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
   {
       bcMAE.addCondition( *bit, condition_type::neumann, &sN );
   }
    gsInfo<<"Boundary conditions:\n"<< bcMAE <<"\n";

    //! [Refinement]
    gsMultiBasis<> dbasis(mpLeft, true);//true: poly-splines (not NURBS)
    
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

    // Set pow for BFO method
    auto IGdim     = G.domainDim();
    // Set factor for BFO method
    auto gammaMAE = factorial(G.domainDim());

    // Set the discretization space
    space u = A.getSpace(dbasis);
    
    // Set the source term with respect to target geometry
    auto ff = A.getCoeff(f, G);

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
    //Initialization of Fast diagonalization solver
    Poisson_FastDiag<double> Poisson(dbasis.basis(0), bcMAE, A.options(), eps);

    u.setup(bcMAE, dirichlet::l2Projection, 0);
    // Compute the system matrix and right-hand side

    // Initialize the system :  identity mapping as initial guess
    A.initSystem();
    // Initialize the system : start Computing the conductivity coeffeicient ...
    // Compute the Neumann terms defined on physical space
    auto g_N = A.getBdrFunction(G);
    auto Neumann_Int{ev.integralBdrBc(bcMAE.get("Neumann"), g_N.tr() * nv(G) )};
    //... nromalisation of density function
    auto CoeffDensity{ev.integral((1.+IntensityMAE*ff.val())* meas(G))};
    auto CoeffConductivity{Neumann_Int/ev.integral(pow(pow(IGdim,IGdim)+gammaMAE * CoeffDensity/(1.+IntensityMAE*ff.val()), 1./IGdim) )};

    setup_time += timer.stop();

    gsInfo<< A.numDofs() <<std::flush;

    timer.restart();
    A.assemble(
    igrad(u, G) * igrad(u, G).tr() * meas(G) + eps * u *u.tr()* meas(G) //matrix
    ,
    u*  CoeffConductivity * (-1.)*pow(pow(IGdim,IGdim)+gammaMAE* CoeffDensity/(1.+IntensityMAE*ff.val()), 1./IGdim) * meas(G) //rhs vector
    );
    
    // Compute the Neumann terms defined on physical space
    //auto g_N = A.getBdrFunction(G);
    A.assembleBdr(bcMAE.get("Neumann"), u * g_N.tr() * nv(G) );
    A.assembleIfc(mp.interfaces(), u.left() * (u_I.tr() * nv(G.left())));
    A.assembleIfc(mp.interfaces(), u.right() * (u_I.tr() * nv(G.right())));

    ma_time += timer.stop();

    gsInfo<< "." <<std::flush;// Assemblying done

    timer.restart();
    // solver.compute( A.matrix() );
    // solVector = solver.solve(A.rhs());
    solVector = Poisson.solve(A.rhs());
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
        //vsolVector = solver.compute(A.matrix()).solve(A.rhs());
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

        A.initSystem(2);
        //Obtain control points for the gradient of mpLeft.comp(Psi)
        A.assemble( v * v.tr() , v * comp.tr() );// blocked by this one
        vsolVector = Poisson.L2ProjectVec(A.rhs());
        //vsolVector = solver.compute(A.matrix()).solve(A.rhs());
        v_sol.extract(PsiLoc);
        //::::::::::::::::::::      end       ::::::::::::::::::::::::: 
        geometryMap PPfLoc = A.getMap(PsiLoc);
        auto ff = A.getCoeff(f, PPfLoc);

        // ...  0  dirichlet for boundaries
        sv0 = solVector;
        u.setup(bcMAE, dirichlet::l2Projection, 0);
    
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
        A.assemble(
        igrad(u, G) * igrad(u, G).tr() * meas(G) +  eps * u * u.tr()* meas(G)//matrix
        ,
        u * CoeffConductivity * (-1.) * ExprMAE * meas(G) //rhs vector
        );

        // Compute the Neumann terms defined on physical space
        auto g_N = A.getBdrFunction(G);
        A.assembleBdr(bcMAE.get("Neumann"), u * g_N.tr() * nv(G) );
        A.assembleIfc(mp.interfaces(), u.left() * (u_I.tr() * nv(G.left())));
        A.assembleIfc(mp.interfaces(), u.right() * (u_I.tr() * nv(G.right())));
        ma_time += timer.stop();

        // gsDebugVar(A.matrix().toDense());
        // gsDebugVar(A.rhs().transpose()   );
        

        gsInfo<< " ." <<std::flush;// Assemblying done

        timer.restart();
        // solver.compute( A.matrix() );
        // solVector = solver.solve(A.rhs());
        solVector = Poisson.solve(A.rhs());

        slv_time += timer.stop();

        gsInfo<< "." <<std::flush; // Linear solving done

        // omp_set_dynamic(0);     // Explicitly disable dynamic teams
        // omp_set_num_threads(1); // Use these threads for later parallel regions

        ++NiterPicard;
        auto l2errRes = math::sqrt(ev.integral( ( igrad(u_lsol,G) - igrad(u_sol,G) ).sqNorm() * meas(G) ));
        if ( l2errRes < tolPicard || ip == maxIter ){
            // ! end Picard loop
            gsInfo<< "\n Niter in Picard : " << NiterPicard << ".. L2 residual : "<<std::scientific<<l2errRes<<"\n";
            break; 
            } // 
    }//for loop
    // omp_set_dynamic(0);     // Explicitly disable dynamic teams
    // omp_set_num_threads(1); // Use these threads for later parallel regions
    gsMultiPatch<> UU;
    u_sol.extract(UU);
    gsWrite(UU, "U_solution");
    auto u_s = A.getCoeff(UU);

    //gsMultiBasis<> gbasis(dbasis);
    //gbasis.reduceContinuity(1);
    space v = A.getSpace(dbasis);
    gsMatrix<> vsolVector;
    solution v_sol = A.getSolution(v, vsolVector);
    A.initSystem(2);

    //gsVector<> pt(2); pt.setConstant(0.5);
    //ev.testEval( v, pt );
    //ev.testEval( igrad(u_sol,G), pt );

    // Obtain control points for the gradient of Psi
    A.assemble( v * v.tr() , v * igrad(u_s,G) );
    vsolVector = solver.compute(A.matrix()).solve(A.rhs());
    gsMultiPatch<> Psi, Psitp;
    v_sol.extract(Psitp);
    //... correct the boundary
    if (PNormalCP)
        ProjectionNormalCPoints(Psitp);

    //::::::::::::::::::::    Compute the composition of geometry maps      :::::::::::::::::::::::::
    // Psi.addAutoBoundaries();
    geometryMap PP = A.getMap(Psitp);
    auto  comp = PP(mpLeft);
    A.initSystem(2);
    //Obtain control points for the gradient of mpLeft.comp(Psi)
    A.assemble( v * v.tr() , v * comp.tr() );// blocked by this one
    vsolVector = solver.compute(A.matrix()).solve(A.rhs());
    v_sol.extract(Psitp);
    Psitp.addAutoBoundaries();
    Psitp.computeTopology();
    gsInfo << "end of adaptive mapping computation\n" << Psitp<< "\n";

    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###   Step 4: Define hierarchical adaptive mapping
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    for(size_t i =0; i<Psitp.nPatches(); ++i)
        Psi.addPatch(gsTHBSpline<2>( dynamic_cast<const gsTensorBSpline<2>&>(Psitp.patch(i)) ));
    Psi.addAutoBoundaries();
    Psi.computeTopology();
    //#-++++++++++++++++++++++++ End of sharing part of any geometry------------------------------
    timer.stop();
    gsInfo<<"\n\nTotal time: "<< setup_time+ma_time+slv_time<<"\n";
    gsInfo<<"     Setup: "<< setup_time <<"\n";
    gsInfo<<"  Assembly: "<< ma_time    <<"\n";
    gsInfo<<"   Solving: "<< slv_time   <<"\n";

    //::::::::::::::::::::   Poisson equation - (manufactured exact solution)         :::::::::::::::::::::::::
    if (true){
    gsBoundaryConditions<> bc;
    bc.setGeoMap(Psi);
    // For simplicity, set Dirichlet boundary conditions
    for ( gsMultiPatch<>::const_biterator
            bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
    {
       bc.addCondition( *bit, condition_type::dirichlet, &s,0, false);
    }
    gsInfo<<"The PDE domain is "<< Psi.detail() << "\n";
    gsInfo<<"Source function is "<< rhs << "\n";
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";


    dbasis.clear();
    gsMultiBasis<> dbasis(Psi, true);//true: poly-splines (not NURBS)

    geometryMap PP = A.getMap(Psi);
    auto ff_GPsi   = A.getCoeff(f, PP);
    
    // --------------- adaptive refinement ---------------
    // Specify cell-marking strategy...
    MarkingStrategy adaptRefCrit = PUCA;
    //MarkingStrategy adaptRefCrit = GARU;
    //MarkingStrategy adaptRefCrit = errorFraction;
    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);
    // Set the discretization space // different boundary condition !
    space ru = A.getSpace(dbasis);

    // Set the source term for Poisson equation
    auto SFunc      = A.getCoeff(rhs, PP);

    // Recover manufactured solution for Poisson equation
    auto u_ex       = ev.getVariable(s, PP);

    // Solution vector and solution variable
    gsMatrix<> rsolVector;
    solution ru_sol = A.getSolution(ru, rsolVector);

    gsVector<>  h1err(numLRefine+1), l2err(numLRefine+1);
    gsVector<int>  DoFPDE(numLRefine+1);
    gsInfo<< "(dot1=assembled, dot2=solved)\n";
    double setup_time(0), ma_time(0), slv_time(0), err_time(0);
    for (int r=0; r<=numLRefine; ++r)
    {
        //::::::::::::::::::::   Poisson equation - (manufactured exact solution)         :::::::::::::::::::::::::
        ru.setup(bc, dirichlet::l2Projection, 0);

        // Compute the system matrix and right-hand side
        // Initialize the system
        A.initSystem();
        setup_time += timer.stop();

        gsInfo<< "Solving PDEs " <<std::flush;
        gsInfo<< A.numDofs() <<std::flush;
        
        //auto h_Tau =  m_h/(2.*coeff_conv.squaredNorm()+m_h);

        timer.restart();
        A.assemble(
        igrad(ru, PP) * igrad(ru, PP).tr() * meas(PP) //matrix
        ,
        ru * SFunc * meas(PP) //rhs vector
        );

        ma_time += timer.stop();

        // gsDebugVar(A.matrix().toDense());
        // gsDebugVar(A.rhs().transpose()   );

        gsInfo<< "." <<std::flush;// Assemblying done

        timer.restart();
        solver.compute( A.matrix() );
        rsolVector = solver.solve(A.rhs());

        slv_time += timer.stop();

        gsInfo<< "." <<std::flush; // Linear solving done

        // omp_set_dynamic(0);     // Explicitly disable dynamic teams
        // omp_set_num_threads(1); // Use these threads for later parallel regions
        DoFPDE[r] = A.numDofs();

        timer.restart();
        l2err[r]= math::sqrt( ev.integral( (u_ex - ru_sol).sqNorm() * meas(PP) ) );

        h1err[r]= math::sqrt(ev.integral( ( igrad(u_ex) - igrad(ru_sol,PP) ).sqNorm() * meas(PP) ));
        err_time += timer.stop();
        gsInfo<< ". " <<std::flush; // Error computations done
        if(r < numLRefine){
        //! [beginRefLoop]
            gsInfo << "====== Loop " << r << " of "
                    <<numLRefine<< " ====adapt Parameter ="<< adaptRefParam << " ======" << "\n";
            // --------------- error estimation/computation ---------------
            // Get the element-wise norms.
            ev.integralElWise( (  ilapl(ru_sol, PP)+ SFunc ).sqNorm() );
            if (IntensityMAE > 1.)
                ev.integralElWise( ff_GPsi );

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
            adaptRefParam = adaptRefParam + FactRefPar;
            }
    }
    //! [Solver loop]    


    timer.stop();
    gsInfo<<"\n\nTotal time: "<< setup_time+ma_time+slv_time+err_time <<"\n";
    gsInfo<<"     Setup: "<< setup_time <<"\n";
    gsInfo<<"  Assembly: "<< ma_time    <<"\n";
    gsInfo<<"   Solving: "<< slv_time   <<"\n";
    gsInfo<<"     Norms: "<< err_time   <<"\n";
    

    //! [Error and convergence rates]
    gsInfo<< "\nDoF_PDE = "<<std::scientific<<DoFPDE.transpose()<<"\n";
    gsInfo<< "L2_error = "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    gsInfo<< "H1_error= "<<std::scientific<<std::setprecision(3)<<h1err.transpose()<<"\n";
    // Assuming DoFPDE, l2err, and h1err are gsMatrix or similar types
    std::ofstream outFile("/home/mbahari/Downloads/error_analysis_ex3.txt", std::ios::app); // Open file in append mode

    if (outFile.is_open())
    {
        outFile << "#DoF_PDE: " << adaptRefParam <<" "<< FactRefPar <<" " << IntensityMAE << " \n"<< std::scientific << DoFPDE.transpose() << "\n";
        outFile << "#L2_error: \n" << std::scientific << std::setprecision(3) << l2err.transpose() << "\n";
        outFile << "#H1_error: \n" << std::scientific << std::setprecision(3) << h1err.transpose() << "\n";
        outFile << "#-------------------------------------------------------------------------------\n"; // Optional separator for readability
        outFile.close(); // Close the file after writing
    }
    else
    {
        gsInfo << "Error: Unable to open file for writing : /home/mbahari/Downloads/error_analysis_ex3.txt.\n";
    }

    //! [Error and convergence rates]
    if (numLRefine>0)
    {
        gsInfo<< "\nEoC (L2): " << std::fixed<<std::setprecision(2)
              <<  ( l2err.head(numRefine).array()  /
                   l2err.tail(numRefine).array() ).log().transpose() / std::log(2.0)
                   <<"\n";

        gsInfo<<   "EoC (H1): "<< std::fixed<<std::setprecision(2)
              <<( h1err.head(numRefine).array() /
                  h1err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";
    }
    //! [Error and convergence rates]
    //! [Export visualization in ParaView]
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        gsParaviewCollection collection("ParaviewOutput/solution", &ev);
        collection.options().setSwitch("plotElements", true);
        collection.options().setSwitch("base64", export_b64);
        collection.options().setInt("plotElements.resolution", 16);
        collection.options().setInt("numPoints", 10000);
        collection.newTimeStep(&Psi);
        collection.addField(ru_sol,"numerical solution");
        collection.addField(igrad(ru_sol,PP),"gradient_numerical solution");
        collection.addField(jac(PP).det(), "Jacobian function");
        collection.addField(u_ex, "exact solution");
        collection.addField(ff_GPsi,"Density function");
        collection.saveTimeStep();
        collection.save();
        gsFileManager::open("ParaviewOutput/solution.pvd");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]
    }

    return EXIT_SUCCESS;


}// end main