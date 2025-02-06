/** @file rh_adaptiveAdvectiondiffusion.cpp

    @brief Tutorial on how to use expression assembler to solve a advection-diffusion equation in rh-adaptive mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris & M. BAHARI
*/

//! [Include namespace]
#include <gismo.h>
# include <gsAssembler/gsAdaptiveRefUtils.h>

using namespace std;
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
    bool plot           = false;
    index_t numRefine   = 3;// for local refinement:  0 means no local h-refinement
    index_t UnifRefine  = 3;// initial refinement: for MAE resolution take at least >=3 for Bejictive mapping 
    index_t DegElevate  = 2; // degree Elevation
    index_t NumArMarEl  = 1; // Number of ring of cells around marked elements
    index_t maxIter     = 50;
    double eps          = 1e-5; // pinalization coefficient
    double tolPicard    = 1e-8;
    double IntensityMAE = 6.;
    real_t adaptRefParam = 0.;     // ... adapt parameter.
    index_t FactRefPar    = 0;    // ... adapt parameter : adaptRefParam += FactRefPar in each iter
    bool ErrorPrint     = true, export_b64 =false;
    // ...PNormalCP: Correct the normal part of the mapping and CornersLshape: adjust the corners of the three patches that form L.
    bool PNormalCP      = true;
    // --------------- adaptive refinement ---------------
    // Specify cell-marking strategy...
    MarkingStrategy adaptRefCrit = PUCA;
    //MarkingStrategy adaptRefCrit = GARU;
    //MarkingStrategy adaptRefCrit = errorFraction;

    gsCmdLine cmd("Tutorial on solving a non-linear Monge-Ampere problem.");
    cmd.addInt("i", "iter", "Maximum number of iterations for the iterative Picard", maxIter);
    cmd.addInt( "e", "DegElevate",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", DegElevate);
    cmd.addInt( "u", "uniformRefine", "Number of Uniform h-refinement loops",  UnifRefine );
    cmd.addInt( "l", "numRefine", "Number of local h-refinement loops",  numRefine );
    cmd.addReal( "a", "adaptRefParam", "parameter for local h-refinement loops",  adaptRefParam );
    cmd.addInt( "p", "FactRefPar", "augement adaptRefParam with such quantity in local h-refinement loops",  FactRefPar );
    cmd.addInt( "c", "NumArMarEl", "augement NumArMarEl with such quantity in local h-refinement loops",  NumArMarEl );
    cmd.addReal( "f", "IntensityMAE", "Intensity of density function",  IntensityMAE);
    cmd.addSwitch( "ErrorPrint", "print Error", ErrorPrint);
    //cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //gsFileData<> fd(fn);
    //gsInfo << "Loaded file "<< fd.lastPath() <<"\n";
    // .... one single patch
    gsMultiPatch<> Psi, mp, mptp = gsNurbsCreator<>::BSplineSquareGrid(1,1,1, 0.0, 0.0);
    // ... need regularity to be at least C^1
    mptp.degreeElevate(DegElevate);
    for(size_t i =0; i<mptp.nPatches(); ++i)
        mp.addPatch(gsTHBSpline<2>( dynamic_cast<const gsTensorBSpline<2>&>(mptp.patch(i)) ));
    mp.addAutoBoundaries();
    
    // // Identity mapping stay fix
    gsFunctionExpr<> sN("x","y",2);

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //..... Test 1 : POISSON EQUATION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    // // Define Stabilization method
    // auto Stabilizationtype = stabilizerCDR::none;
    // // convection coefficient
    // gsFunctionExpr<> coeff_conv("0.","0.",2);
    // // diffusion coefficient:
    // gsFunctionExpr<> coeff_diff("1.","0","0","1.",2);
    // // For a posterior error estimate
    // gsFunctionExpr<> coeff_diffMax("1.",2);
    // // reaction coefficient:
    // gsFunctionExpr<> coeff_reac("0",2);
    /* *********************** Test 1 *********************** */
    // // Define  Dirichlet boundary conditions
    // gsFunctionExpr<> Dg("1./(1.+exp((y - x  - 0.2)/0.01))", 2);
    // // Manufactured solition
    // gsFunctionExpr<> s("1./(1.+exp((y - x  - 0.2)/0.01))",2);
    // // // Right-hand side function
    // gsFunctionExpr<> SourceFunc("4.12230724487712e-5*exp(-100.0*x + 100.0*y)/(2.06115362243856e-9*exp(-100.0*x + 100.0*y) + 1.0)**2 - 1.69934170211664e-13*exp(-200.0*x + 200.0*y)/(2.06115362243856e-9*exp(-100.0*x + 100.0*y) + 1.0)**3",2);
    // // Manufactured density function
    // gsFunctionExpr<> f("( 1/(1.+exp((y -x  - 0.3)/0.01)) - 1/(1.+exp((y - x  - 0.1)/0.01))  )",2);
    /* *********************** Test 2 *********************** */
    // // Define  Dirichlet boundary conditions
    // gsFunctionExpr<> Dg("1./(1.+exp(100 * ( x**2 + (y-0.5)**2-0.75*sin(pi*y)) ))", 2);
    // // Manufactured solition
    // gsFunctionExpr<> s("1./(1.+exp(100 * ( x**2 + (y-0.5)**2-0.75*sin(pi*y)) ))",2);
    // // // Right-hand side function
    // gsFunctionExpr<> SourceFunc("40000.0*x**2*exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y))/(exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y)) + 1.0)**2 - 80000.0*x**2*exp(200*x**2 + 200*(y - 0.5)**2 - 150.0*sin(pi*y))/(exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y)) + 1.0)**3 + 1.0*(75.0*pi**2*sin(pi*y) + 200)*exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y))/(exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y)) + 1.0)**2 + 1.0*(40000*(y - 0.375*pi*cos(pi*y) - 0.5)**2)*exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y))/(exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y)) + 1.0)**2 + 200.0*exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y))/(exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y)) + 1.0)**2 - 80000.0*(y - 0.375*pi*cos(pi*y) - 0.5)**2*exp(200*x**2 + 200*(y - 0.5)**2 - 150.0*sin(pi*y))/(exp(100*x**2 + 100*(y - 0.5)**2 - 75.0*sin(pi*y)) + 1.0)**3",2);
    // // Manufactured density function
    // gsFunctionExpr<> f("exp(-90 * ( x**2 + (y-0.5)**2-0.75*sin(pi*y))**2 )",2);

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // //..... Test 2 ADVECTION DUFFFUSION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    // Define Stabilization method
    auto Stabilizationtype = stabilizerCDR::SUPG;
    // Define  Dirichlet boundary conditions
    gsFunctionExpr<> Dg("if( y <= 0.2*(1.-x), 1,0)", 2);
    // Manufactured solition
    gsFunctionExpr<> s("if( y <1., 1./(1.+exp((y - x  - 0.2)/0.0001))-1/(1.+exp((1.-x)/0.0001)),0)",2);
    // convection coefficient:
    gsFunctionExpr<> coeff_conv("cos(pi/4)","sin(pi/4)",2);
    // diffusion coefficient:
    gsFunctionExpr<> coeff_diff("0.000001","0","0","0.000001",2);
    // For a posterior error estimate
    gsFunctionExpr<> coeff_diffMax("0.000001",2);
    // reaction coefficient:
    gsFunctionExpr<> coeff_reac("0",2);
    // // Right-hand side function
    gsFunctionExpr<> SourceFunc("0.",2);
    //Manufactured density function 1./cosh(100. * ( -x - 0.2 + y ))
    gsFunctionExpr<> f("( 1./cosh( 10.*( -x+y -0.2 ) )**2 + 1/(1.+exp((0.95-x)/0.01)) )",2);

    gsInfo<<"Source function "<< f << "\n";

    gsInfo<<"The Initial domain is "<< mp.detail() << "\n";

    gsBoundaryConditions<> bc_mae;
    bc_mae.setGeoMap(mp);
    // For simplicity, set Neumann boundary conditions
    for ( gsMultiPatch<>::const_biterator
            bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
    {
       bc_mae.addCondition( *bit, condition_type::neumann, &sN,0, false);
    }
    //gsDebugVar( bc.allConditions()[0].parametric() );
    gsInfo<<"Boundary conditions:\n"<< bc_mae <<"\n";

    //gsOptionList Aopt;

    //! [Refinement]
    gsMultiBasis<> dbasis(mp, true);//true: poly-splines (not NURBS)

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    //dbasis.setDegree( dbasis.maxCwiseDegree() + numElevate);
    //gsInfo << dbasis.degree(0) << " degree  \n";

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
#ifdef _OPENMP
    gsInfo<< "Available threads: "<< omp_get_max_threads() <<"\n";
#endif

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    //! [Solver loop]
    gsSparseSolver<>::CGDiagonal solver;

    gsVector<>  h1err(numRefine+1), l2err(numRefine+1);
    gsVector<int>  DoFPDE(numRefine+1);
    gsInfo<< "(dot1=assembled, dot2=solved, dot3=nonlinear_loop,dot4=got_error)\n";
    double setup_time(0), ma_time(0), slv_time(0), err_time(0);    
    gsStopwatch timer;

     //! [Refinement]
    mp.uniformRefine();
    for (int r=0; r<=UnifRefine; ++r)
    {
        dbasis.uniformRefine();
        //mp.uniformRefine();
        Psi.uniformRefine();
    }
    //::::::::::::::::::::      mesh adaptation solver         :::::::::::::::::::::::::

    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    //A.setOptions(Aopt);
    gsInfo<<"Active options:\n"<< A.options() <<"\n";

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);
    auto tt = A.rhs();
    // Set the discretization space
    space u = A.getSpace(dbasis); 

    // Set the geometry map
    geometryMap G = A.getMap(mp);

    // Set pow for BFO method
    auto IGdim     = G.domainDim();
    // Set factor for BFO method
    auto gammaMAE = factorial(G.domainDim());

    // Set the source term
    auto ff = A.getCoeff(f, G);

    // Solution vector and solution variable
    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);
    
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###   Step 2: r_refinement
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    Poisson_FastDiag<double> Poisson(dbasis.basis(0), bc_mae, A.options(), eps);

    //... nromalisation of density function
    auto CoeffDensity{ev.integral((1.+IntensityMAE*ff.val()) * meas(G))};
    // Initialize the system : start Computing the conductivity coeffeicient ...
    // Compute the Neumann terms defined on physical space
    auto g_N = A.getBdrFunction(G);
    auto Neumann_Int{ev.integralBdrBc(bc_mae.get("Neumann"), g_N.tr() * nv(G) )};
    // ...
    auto CoeffConductivity{Neumann_Int/ev.integral(pow(IGdim*IGdim+gammaMAE * CoeffDensity/(1.+IntensityMAE*ff.val()), 1./IGdim) * meas(G))};
    //... end 

    //u.setup(bc_mae, dirichlet::l2Projection, 0);
    // Initialize the system :  identity mapping as initial guess
    A.initSystem();
    setup_time += timer.stop();

    gsInfo<< "\nDoFs: " << A.numDofs() <<std::flush << "\n";

    timer.restart();

    A.assemble(
    igrad(u, G) * igrad(u, G).tr() * meas(G) + eps * u *u.tr()* meas(G) //matrix
    ,
    u*  CoeffConductivity * (-1.)*pow(pow(IGdim,IGdim)+gammaMAE * CoeffDensity/(1.+IntensityMAE*ff.val()), 1./IGdim) * meas(G) //rhs vector
    );
    
    // Compute the Neumann terms defined on physical space
    A.assembleBdr(bc_mae.get("Neumann"), u * g_N.tr() * nv(G) );

    ma_time += timer.stop();

    // gsDebugVar(A.matrix().toDense());
    // gsDebugVar(A.rhs().transpose()   );

    gsInfo<< "." <<std::flush;// Assemblying done

    timer.restart();
    solVector = Poisson.solve(A.rhs());
    // solver.compute( A.matrix() );
    // solVector = solver.solve(A.rhs());
    slv_time += timer.stop();

    gsInfo<< "." <<std::flush; // Linear solving done

    // Picard loop
    gsMatrix<> sv0; //
    solution u_lsol = A.getSolution(u, sv0);
    for(int ip{0}; ip<=maxIter; ++ip)
    {
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
        // vsolVector = solver.compute(A.matrix()).solve(A.rhs());
        vsolVector = Poisson.L2ProjectVec(A.rhs());

        v_sol.extract(Psi);                
        
        // Set the geometry optimal map
        geometryMap PP = A.getMap(Psi);
        auto ff = A.getCoeff(f,PP);

        // ...  0  dirichlet for boundaries
        sv0 = solVector;
        
        //u.setup(bc_mae, dirichlet::l2Projection, 0);

        // Initialize the system
        A.initSystem();
        setup_time += timer.stop();
        timer.restart();
        // Compute the system matrix and right-hand side ... Monge-Ampere eqaution .....
        auto  ExprMAE     = pow( pow(div(PP).val(),IGdim) + gammaMAE*(CoeffDensity/(1.+IntensityMAE*ff.val()) - jac(PP).det()), 1./IGdim);
        auto IntegDensity = ev.integral(ExprMAE);
        
        // .. update Coeffeicient of conductivity
        CoeffConductivity = Neumann_Int/IntegDensity;
        // MAE system
        A.assemble(
        igrad(u, G) * igrad(u, G).tr() * meas(G) +  eps * u * u.tr()* meas(G) //matrix
        ,
        u * CoeffConductivity * (-1.) * ExprMAE * meas(G) //rhs vector
        );

        // Compute the Neumann terms defined on physical space
        auto g_N = A.getBdrFunction(G);
        A.assembleBdr(bc_mae.get("Neumann"), u * g_N.tr() * nv(G) );

        ma_time += timer.stop();            

        gsInfo<< " ." <<std::flush;// Assemblying done

        timer.restart();
        solVector = Poisson.solve(A.rhs());
        // solver.compute( A.matrix() );
        // solVector = solver.solve(A.rhs());
        
        slv_time += timer.stop();

        gsInfo<< "." <<std::flush; // Linear solving done

        // omp_set_dynamic(0);     // Explicitly disable dynamic teams
        // omp_set_num_threads(1); // Use these threads for later parallel regions

        auto l2errRes = math::sqrt(ev.integral( ( igrad(u_lsol,G) - igrad(u_sol,G) ).sqNorm() * meas(G) ));
        if ( l2errRes < tolPicard || ip == maxIter ){
            // ! end Picard loop
            gsInfo<< "\n Niter in Picard : " << ip << ".. L2 residual : "<<std::scientific<<l2errRes<<"\n";
            break; 
            } // 
    }//for loop
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###   Step 3: Correct boundary
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    if (PNormalCP)
        ProjectionNormalCPoints(Psi);
    //  construct boundaries and interfaces
    if(mp.nPatches()>1){
    Psi.addInterface(0,2,1,1);
    Psi.addInterface(0,4,2,3);
    Psi.addAutoBoundaries();
    //Psi.computeTopology();
    }
    Psi.addAutoBoundaries();

    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###   Step 4: Define hierarchical adaptive mapping
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    //gsMultiPatch<> Psi = Psi;
    // for(size_t i =0; i<Psi.nPatches(); ++i)
    //     Psi.addPatch(gsTHBSpline<2>( dynamic_cast<const gsTensorBSpline<2>&>(Psi.patch(i)) ));
    // //Psi.addAutoBoundaries();
    // Psi.computeTopology();
    gsInfo<<"The PDE domain is "<< Psi.detail() << "\n";

    // Set Dirichlet boundary conditions
    gsBoundaryConditions<> bc;
    bc.setGeoMap(Psi);
    // given by exact solution Dg on all boundaries:
    for ( gsMultiPatch<>::const_biterator
                bit = Psi.bBegin(); bit != Psi.bEnd(); ++bit)
    {
        bc.addCondition( *bit, condition_type::dirichlet, &Dg );
    }
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";

    //! [Refinement]
    //gsMultiBasis<> hdbasis(Psi, true);//true: poly-splines (not NURBS)

    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###   Step 5: local h-refinement
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    gsInfo << "Patches: "<< Psi.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";

    // // Set the geometry optimal map
    geometryMap PP    = ev.getMap(Psi);  

    // Recover manufactured solution for Poisson equation
    auto u_ex         = ev.getVariable(s, PP);
    auto rho          = ev.getVariable(f, PP);

    // Recover rhs for Poisson equation
    auto SFunc        = ev.getVariable(SourceFunc, PP);
    // Coeffs for advection-reaction diffusion equation
    auto coeff_convPP = ev.getVariable(coeff_conv, PP);
    auto coeff_diffPP = ev.getVariable(coeff_diffMax, PP);
    auto coeff_reacPP = ev.getVariable(coeff_reac, PP);
    // numerical solutionas Vector
    gsMatrix<> rsolVector;
    for (int r=0; r<=numRefine; ++r)
    {
        // --------------- define Pde ---------------
        timer.restart();
        //! [definePde]
        gsConvDiffRePde<real_t> cdrPde(Psi, bc, & coeff_diff,& coeff_conv, & coeff_reac, & SourceFunc);
        //! [definePde]
        //! [constructAssembler]
        // Construct assembler
        gsCDRAssembler<real_t> cdrAss( cdrPde, dbasis);
        // Set stabilization flag to 1 = SUPG
        cdrAss.options().setInt("Stabilization", Stabilizationtype);
        // Compute Dirichlet values by L2-projection
        // Caution: Interpolation does not work for locally refined (T)HB-splines!
        cdrAss.options().setInt("DirichletValues",dirichlet::l2Projection);
        //! [constructAssembler]
        gsInfo<< "." <<std::flush;// Assemblying done
        ma_time += timer.stop();
        // Generate system matrix and load vector
        cdrAss.assemble();
        // Solution vector and solution variable
        timer.restart();
        // Solve the system
        rsolVector = gsSparseSolver<>::BiCGSTABILUT( cdrAss.matrix() ).solve( cdrAss.rhs() );

        slv_time += timer.stop();
        gsInfo<< "DoFs in PDEs " << cdrAss.numDofs() <<std::flush;
        DoFPDE[r] = cdrAss.numDofs();
        gsField<> solField;        // Construct the solution as a scalar field
        solField = cdrAss.constructSolution(rsolVector);

        ev.setIntegrationElements(dbasis);

        gsExprEvaluator<>::variable ru_sol = ev.getVariable(solField.fields());

        // omp_set_dynamic(0);     // Explicitly disable dynamic teams
        // omp_set_num_threads(1); // Use these threads for later parallel regions

        //! [errorComputation]
        timer.restart();
        l2err[r]= math::sqrt( ev.integral( (u_ex - ru_sol).sqNorm() * meas(PP) ) );
        h1err[r]= math::sqrt(ev.integral( ( igrad(u_ex) - igrad(ru_sol,PP) ).sqNorm() * meas(PP) ));
        err_time += timer.stop();
        gsInfo<< ". " <<std::flush; // Error computations done

        if(r < numRefine){
        //! [beginRefLoop]
            gsInfo << "====== Loop " << r << " of "
                    <<numRefine<< " ====adapt Parameter ="<< adaptRefParam << " ======" << "\n";
        // --------------- error estimation/computation ---------------
        // for test :ev.integralElWise( ( ilapl(ru_sol, PP) + SFunc ).sqNorm()*meas(PP) );
        // Get the element-wise norms.
        ev.integralElWise( ( coeff_diffPP * ilapl(ru_sol,PP) - igrad(ru_sol, PP)[0]*0.7071067811865476- igrad(ru_sol, PP)[1]*0.7071067811865476 - coeff_reacPP * ru_sol + SFunc).sqNorm() );
        if(IntensityMAE > 1.)
            ev.integralElWise( rho );

        const std::vector<real_t> eltErrs  = ev.elementwise();
        //! [errorComputation]

        //! [adaptRefinementPart]
        // Mark elements for refinement, based on the computed local errors and
        // the refinement-criterion and -parameter.
        std::vector<bool> elMarked( eltErrs.size() );
        gsMarkElementsForRef( eltErrs, adaptRefCrit, adaptRefParam, elMarked);
        gsInfo <<"Marked "<< std::count(elMarked.begin(), elMarked.end(), true) <<" elements.\n";
        // Refine the marked elements with a 1-ring of cells around marked elements
        gsRefineMarkedElements( dbasis, elMarked, NumArMarEl);
        gsRefineMarkedElements( Psi, elMarked, NumArMarEl);

       // Call repair interfaces to make sure that the new meshes
       // match along patch interfaces.
       dbasis.repairInterfaces( Psi.interfaces() );
       //! [refreshAssembler]
       cdrAss.refresh();
       if(r%2==0)
       NumArMarEl = NumArMarEl + FactRefPar;
        }
    if(plot && r == numRefine){
        // gsInfo<<"Storing paraview...\n";
        // // Write the computed solution to paraview files
        // gsWriteParaview<>( solField, "adaptRef", 1000, true);
        gsInfo<<"Making in Paraview...\n";
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
        collection.addField(rho,"Density function");
        collection.saveTimeStep();
        collection.save();
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

    if (ErrorPrint && numRefine>0)
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
       //! [Export to Paraview]
       // Export the final solution
    if(plot){
        //------------------------------------
        gsInfo<<"Plotting in Paraview...\n";
        // Run paraview
        gsFileManager::open("ParaviewOutput/solution.pvd");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;


}// end main