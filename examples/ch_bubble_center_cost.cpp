/** @file ch_circle_time_linear_load.cpp

    @brief Tutorial on how to use expression assembler to solve the Cahn-Hilliard equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Marsala (UniFi)
               H.M. Verhelst (UniFi)
               L. Venta Viñuela (UniPv)


    Tensor product (check flags in .xml file): with lambda = 0.1
    //  =========== QUADRATIC BASIS ===========
    ./bin/ch_circle_time_linear_load --plot -r 6 -t 0.001 -N 20 -c 0 -v 2 -s 0 -f pde/cahn_hilliard_bvp.xml -o "hola"
    ./bin/ch_circle_time_linear_load --plot -r 7 -t 0.001 -N 20 -c 0 -v 2 -s 0 -f pde/cahn_hilliard_bvp.xml -o "hola"
    ./bin/ch_circle_time_linear_load --plot -r 8 -t 0.001 -N 20 -c 0 -v 2 -s 0 -f pde/cahn_hilliard_bvp.xml -o "hola"

    and clamped BC for the essential flux boundary condition
    -----------------------------------------------------------------------
    TODO;
    - Change hmax to a gsExprAssembler<>::element el; el.diam();
    -----------------------------------------------------------------------

*/

//! [Include namespace]
#include <gismo.h>

using namespace gismo;
//! [Include namespace]

template <class T>
T interpolate_crs(const gsMultiBasis<T> & dbasis_fine,
                  const gsMultiBasis<T> & dbasis_coarse,
                  const gsGeometry<T> & fine,
                  gsMatrix<T> & coefs,
                  gsMultiPatch<T> & mp,
                  index_t & projection_Crs)
{
    if (projection_Crs == 0)
    {
        gsL2Projection<real_t>::project(dbasis_fine, dbasis_coarse,fine,mp,coefs);
    }
    else if(projection_Crs == 1)
    {
        gsQuasiInterpolate<real_t>::localIntpl(dbasis_coarse.basis(0),fine,coefs);
    }

    gsExprEvaluator<> ev;
    auto G = ev.getMap(mp);
    ev.setIntegrationElements(dbasis_fine);

    gsGeometry<>::uPtr coarse_ptr = dbasis_coarse.basis(0).makeGeometry(coefs);
    auto cfine = ev.getVariable(fine,G);
    auto ccoarse = ev.getVariable(*coarse_ptr,G);

    real_t err_proj_crs = ev.integral(meas(G) * (ccoarse-cfine).sqNorm());

    return err_proj_crs;
}

template <class T>
T interpolate_ref(const gsMultiBasis<T> & dbasis_fine,
                  const gsGeometry<T> & coarse,
                  gsMatrix<T> & coefs,
                  gsMultiPatch<T> & mp)
{
    gsQuasiInterpolate<real_t>::localIntpl(dbasis_fine.basis(0),coarse,coefs);

    gsExprEvaluator<> ev;
    auto G = ev.getMap(mp);
    ev.setIntegrationElements(dbasis_fine);

    gsGeometry<>::uPtr fine_ptr = dbasis_fine.basis(0).makeGeometry(coefs);
    auto cfine = ev.getVariable(*fine_ptr,G);
    auto ccoarse = ev.getVariable(coarse,G);

    real_t err_proj_ref = ev.integral(meas(G) * (ccoarse-cfine).sqNorm());

    return err_proj_ref;
}

template <short_t dim, class T>
void solve( gsMultiPatch<T> & mp,
            gsFunctionExpr<T> & source,
            gsBoundaryConditions<T> & bc,
            gsOptionList & CHopt,
            gsOptionList & TIMEopt,
            gsOptionList & MESHopt,
            gsOptionList & Aopt,
            real_t & dt,
            index_t & maxSteps,
            index_t & plotmod,
            bool & plot,
            bool & plot_error,
            index_t & numRefine,
            index_t & numElevate,
            index_t & verbose,
            bool & random,
            index_t & projection_Crs,
            index_t & pattern,
            std::string out)
{

    real_t assemblyTime = 0;
    real_t solverTime = 0;
    real_t projectionTime = 0;
    index_t nSolves = 0;
    index_t nSolvesStep = 0;

    gsStopwatch clock;

    real_t theta    = CHopt.askReal("theta",1.5);
    real_t lambda   = CHopt.askReal("lambda",1/(32*pow(EIGEN_PI,2)));
    real_t M0       = CHopt.askReal("M0",0.005);
    real_t penalty  = 1e4*lambda;

    //! [Prepare the basis]
    gsMultiBasis<> dbasis_tmp(mp,true);
    gsMultiBasis<> dbasis;

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis_tmp.setDegree( dbasis_tmp.maxCwiseDegree() + numElevate);

    // gsFileData<> fd1;
    // fd1.read("basis_cir_mark_"+std::to_string(MESHopt.getReal("CoarsenParam"))+"_r_"+std::to_string(numRefine)+".xml");
    // fd1.getId(0,dbasis);
    // gsInfo << "Loaded "<< "basis_cir_mark_"+std::to_string(MESHopt.getReal("CoarsenParam"))+"_r_"+std::to_string(numRefine)+".xml" <<"\n";

    if (MESHopt.askSwitch("Adaptive",true))         // Load the basis from file
    {     
        gsFileData<> fd1;
        fd1.read("basis_circle_mark_"+std::to_string(MESHopt.getReal("CoarsenParam"))+"_r_"+std::to_string(numRefine)+".xml");
        fd1.getId(0,dbasis);
        gsInfo << "Loaded "<<  "basis_circle_mark_"+std::to_string(MESHopt.getReal("CoarsenParam"))+"_r_"+std::to_string(numRefine)+".xml" <<"\n";
        MESHopt.setSwitch("THB",true);
    }
    else // Tensor product basis
    {
        if (MESHopt.askSwitch("THB",true))
        {
            // Cast every basis of dbasis to a gsTHBSplineBasis
            for (size_t p=0; p!=dbasis_tmp.nBases(); p++)
            {
                // TODO: Make dimension-independent over the template
                if (gsTensorBSplineBasis<dim,real_t> * b = dynamic_cast<gsTensorBSplineBasis<dim,real_t>*>(&dbasis_tmp.basis(p)))
                    dbasis.addBasis(new gsTHBSplineBasis<dim,real_t>(*b));
                else if (gsTHBSplineBasis<dim,real_t> * b = dynamic_cast<gsTHBSplineBasis<dim,real_t>*>(&dbasis_tmp.basis(p)))
                    dbasis.addBasis(b->clone());
                else
                    GISMO_ERROR("Basis is neither a gsTHBSplineBasis nor a gsTensorBSplineBasis");
    
                // Refine the basis for `numRefine` levels
                gsMatrix<> box = dbasis.basis(p).support();
                for (index_t r = 0; r!=numRefine; r++)
                    dbasis.basis(p).refine(box);
            }
        }
        else
        {
            for (size_t p=0; p!=dbasis_tmp.nBases(); p++)
            {
                dbasis.addBasis(dbasis_tmp.basis(p).clone());
                for (index_t r = 0; r!=numRefine; r++)
                    dbasis.basis(p).uniformRefine();
            }
        }    
    }

    

    // Determine maximum mesh size
    real_t hmax = 0;
    for (size_t p=0; p!=dbasis.nBases(); p++)
    {
        hmax = math::max(hmax, dbasis.basis(p).getMaxCellLength());
        gsInfo<<"Maximum degree of the basis: "<<dbasis.basis(p).maxDegree()<<"\n";
    }

    if (random)
        if(pattern==0)
            out = out + "_random_nuc_N_" + std::to_string(maxSteps) + "_dt_" + std::to_string(dt) + "_lambda_" + std::to_string(lambda) + "_r_" + std::to_string(numRefine) + "_degree_" + std::to_string(dbasis.basis(0).maxDegree()) + "_prjCrs_" + std::to_string(projection_Crs) + "_THB_" + std::to_string(MESHopt.askSwitch("THB")) + "_Adaptive_" + std::to_string(MESHopt.askSwitch("Adaptive")) + "_refIt_" + std::to_string(MESHopt.askInt("RefIt"));
        else
            out = out + "_random_spin_N_" + std::to_string(maxSteps) + "_dt_" + std::to_string(dt) + "_lambda_" + std::to_string(lambda) + "_r_" + std::to_string(numRefine) + "_degree_" + std::to_string(dbasis.basis(0).maxDegree()) + "_prjCrs_" + std::to_string(projection_Crs) + "_THB_" + std::to_string(MESHopt.askSwitch("THB")) + "_Adaptive_" + std::to_string(MESHopt.askSwitch("Adaptive"))+ "_refIt_" + std::to_string(MESHopt.askInt("RefIt"));
    else
        out = out + "_analytical_N_" + std::to_string(maxSteps) + "_dt_" + std::to_string(dt) + "_lambda_" + std::to_string(lambda) + "_r_" + std::to_string(numRefine) + "_degree_" + std::to_string(dbasis.basis(0).maxDegree()) + "_prjCrs_" + std::to_string(projection_Crs) + "_THB_" + std::to_string(MESHopt.askSwitch("THB")) + "_Adaptive_" + std::to_string(MESHopt.askSwitch("Adaptive")) + "_refIt_" + std::to_string(MESHopt.askInt("RefIt"));

    //! [Prepare the basis]

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // Reals for time integration

    // Generalized-alpha method parameters
    real_t rho_inf = TIMEopt.askReal("rho_inf",0.5);
    real_t alpha_m = 0.5*(3-rho_inf) / (1+rho_inf);
    real_t alpha_f = 1 / (1+rho_inf);
    real_t gamma   = 0.5 + alpha_m - alpha_f;
    real_t dt_old = dt;
    real_t t_rho = TIMEopt.askReal("t_rho",0.9);
    real_t t_err = 1;
    index_t lmax = 1;

    real_t tmp_alpha_m = 1;
    real_t tmp_alpha_f = 1;
    real_t tmp_gamma   = 1;
    // time stepping options
    index_t maxIt = 50;



    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    A.options().setSwitch("SameElement",Aopt.askSwitch("SameElement",true));

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    // A.setIntegrationElements(dbasis);

    // ========= NEW ADD LINEAR BASIS FOR L2 PROJECTION OF dc =========
    // gsMultiBasis<> dbasis_tmp_2(mp,true);
    // dbasis_tmp_2.setDegree(3);
    // gsMultiBasis<> dbasis_lin;
    // if (MESHopt.askSwitch("THB",true))
    // {
    //     // Cast every basis of dbasis to a gsTHBSplineBasis
    //     for (size_t p=0; p!=dbasis_tmp_2.nBases(); p++)
    //     {
    //         // TODO: Make dimension-independent over the template
    //         if (gsTensorBSplineBasis<dim,real_t> * b = dynamic_cast<gsTensorBSplineBasis<dim,real_t>*>(&dbasis_tmp_2.basis(p)))
    //             dbasis_lin.addBasis(new gsTHBSplineBasis<dim,real_t>(*b));
    //         else if (gsTHBSplineBasis<dim,real_t> * b = dynamic_cast<gsTHBSplineBasis<dim,real_t>*>(&dbasis_tmp_2.basis(p)))
    //             dbasis_lin.addBasis(b->clone());
    //         else
    //             GISMO_ERROR("Basis is neither a gsTHBSplineBasis nor a gsTensorBSplineBasis");

    //         // Refine the basis for `numRefine` levels
    //         gsMatrix<> box = dbasis_lin.basis(p).support();
    //         for (index_t r = 0; r!=numRefine; r++)
    //             dbasis_lin.basis(p).refine(box);
    //     }
    // }

    // // Refine the basis for `numRefine` levels
    // gsMatrix<> box = dbasis_lin.basis(0).support();
    // for (index_t r = 0; r!=1; r++)
    //     dbasis_lin.basis(0).refine(box);

    // gsDebugVar(dbasis_lin.basis(0));
    // gsDebugVar(dbasis_lin.basis(0).maxDegree());
    
    // =====================================================================

    gsMultiBasis<> ibasis;
    // Cast every basis of dbasis to a gsTHBSplineBasis
    for (size_t p=0; p!=dbasis.nBases(); p++)
    {
        GISMO_ASSERT((dynamic_cast<gsTHBSplineBasis<dim,real_t>*>(&dbasis.basis(p))), "Basis is not a gsTHBSplineBasis");
        gsTHBSplineBasis<dim,real_t> & b = static_cast<gsTHBSplineBasis<dim,real_t> & >(dbasis.basis(p));
        ibasis.addBasis(b.tensorLevel(b.maxLevel()).clone());
    }
    A.setIntegrationElements(ibasis); // USE FINEST TENSOR BASIS FOR INTEGRATION
    gsExprEvaluator<> ev(A);
    A.options().update(A.options(),gsOptionList::addIfUnknown); 
    ev.setIntegrationElements(dbasis);

    gsDebugVar(dbasis.basis(0).maxDegree());
    gsDebugVar(ibasis.basis(0));

    // Set the geometry map
    geometryMap G = A.getMap(mp);

    // Set the discretization space
    space w = A.getSpace(dbasis);

    // basis.init(dbasis, cf);

    // Solution vector and solution variable
    gsMatrix<> Cnew, Calpha, Cold, Cold1;
    gsMatrix<> dCnew,dCalpha,dCold, dCupdate;
    // Solution vectors with eliminated and free DoFs
    gsMatrix<> CnewF, dCnewF;

    // Solution variables for the previous and next solutions (before and after time step)
    gsMultiPatch<> mp_cold, mp_dcold, mp_cnew, mp_dcnew;
    auto     cold  = A.getCoeff(mp_cold );
    auto     dcold = A.getCoeff(mp_dcold);
    auto     cnew  = A.getCoeff(mp_cnew );
    auto     dcnew = A.getCoeff(mp_dcnew);
    // New solution variables for source term;
    gsMultiPatch<> mp_qnew, mp_qold;
    auto     qnew  = A.getCoeff(mp_qnew );
    auto     qold  = A.getCoeff(mp_qold );

    solution cnew_sol = A.getSolution(w, Cnew); // Cnew
    solution dcnew_sol = A.getSolution(w, dCnew); // \dot{Cnew}

    // Solution variables for the intermediate solutions (during time integration)
    gsConstantFunction<> tmp_alpha_f_func(tmp_alpha_f,dim);
    gsConstantFunction<> tmp_alpha_m_func(tmp_alpha_m,dim);
    auto     af = A.getCoeff(tmp_alpha_f_func);
    auto     am = A.getCoeff(tmp_alpha_m_func);

    auto     c  = cold  + af * ( cnew_sol -  cold);
    auto     dc = dcold + am * (dcnew_sol - dcold);
    auto  gradc = igrad(cold,G) + af * (igrad(cnew_sol,G) - igrad(cold,G));
    auto  laplc = ilapl(cold,G) + af * (ilapl(cnew_sol,G) - ilapl(cold,G));

    // Source term
    auto     qq  = qold  + af * (qnew -  qold);

    // Source (volume integral) function for manufactured solution cos(6*pi*x) * cos(6*pi*y) * cos(2/3*pi*t)
    // gsFunctionExpr<> sourceQold("(20*(tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18) + 1)*(192*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18) - 192)*(((4*(3*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)^2 - 1)*(2*x^2 - 2*x + 2*y^2 - 2*y + 1) - (3*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)^2 + 120*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)*(3*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)^2 - 2)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) - 1)*(x^2/2 - x/2 + y^2/2 - y/2 + 1/4) - (tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2))/30 + 240*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)*(3*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)^2 - 2)*(2*x^2 - 2*x + 2*y^2 - 2*y + 1)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 1/1800)*(x^2 - x + y^2 - y + 1/2)^2)/2 + (tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2))/750000000 - ((15*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)^2 + 120*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)*(3*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)^2 - 2)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) - 5)*(x^2/4 - x/4 + y^2/4 - y/4 + 1/8))/100000000 - (3*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)^2)/2000000000000 + (3*(x^2/3 - x/3 + y^2/3 - y/3 + 1/6)*(4*(3*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)^2 - 1)*(5*x^2 - 5*x + 5*y^2 - 5*y + 5/2) - (3*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)^2 + 60*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)*(3*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)^2 - 2)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) - 1)*(x^2 - x + y^2 - y + 1/2) - (2*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2))/15 + 120*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)*(3*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)^2 - 2)*(2*x^2 - 2*x + 2*y^2 - 2*y + 1)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 1/450))/10000 + ((3*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)^2 - 1)*(2*x^2 - 2*x + 2*y^2 - 2*y + 1))/12500000 + (pi*sin((pi*t)/3)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2)*(x^2 - x + y^2 - y + 5001/10000)^3)/480 - 391/18000000000000))/((4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2)*(2*x^2 - 2*x + 2*y^2 - 2*y + 5001/5000)^3)",2);
    // gsFunctionExpr<> sourceQnew("(20*(tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18) + 1)*(192*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18) - 192)*(((4*(3*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)^2 - 1)*(2*x^2 - 2*x + 2*y^2 - 2*y + 1) - (3*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)^2 + 120*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)*(3*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)^2 - 2)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) - 1)*(x^2/2 - x/2 + y^2/2 - y/2 + 1/4) - (tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2))/30 + 240*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)*(3*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)^2 - 2)*(2*x^2 - 2*x + 2*y^2 - 2*y + 1)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 1/1800)*(x^2 - x + y^2 - y + 1/2)^2)/2 + (tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2))/750000000 - ((15*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)^2 + 120*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)*(3*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)^2 - 2)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) - 5)*(x^2/4 - x/4 + y^2/4 - y/4 + 1/8))/100000000 - (3*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)^2)/2000000000000 + (3*(x^2/3 - x/3 + y^2/3 - y/3 + 1/6)*(4*(3*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)^2 - 1)*(5*x^2 - 5*x + 5*y^2 - 5*y + 5/2) - (3*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)^2 + 60*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)*(3*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)^2 - 2)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) - 1)*(x^2 - x + y^2 - y + 1/2) - (2*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2))/15 + 120*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)*(3*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)^2 - 2)*(2*x^2 - 2*x + 2*y^2 - 2*y + 1)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 1/450))/10000 + ((3*tanh(3*cos((pi*t)/3) - 30*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 18)^2 - 1)*(2*x^2 - 2*x + 2*y^2 - 2*y + 1))/12500000 + (pi*sin((pi*t)/3)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2)*(x^2 - x + y^2 - y + 5001/10000)^3)/480 - 391/18000000000000))/((4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2)*(2*x^2 - 2*x + 2*y^2 - 2*y + 5001/5000)^3)",2);
    // gsFunctionExpr<> sourceQold("(40*(192*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12) - 192)*(tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12) + 1)*((3*(4*(3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 1)*(5*x^2 - 5*x + 5*y^2 - 5*y + 5/2) - (3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 + 40*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 2)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) - 1)*(x^2 - x + y^2 - y + 1/2) - (tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2))/5 + 80*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 2)*(2*x^2 - 2*x + 2*y^2 - 2*y + 1)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 1/200)*(x^2/3 - x/3 + y^2/3 - y/3 + 1/6))/10000 - ((15*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 + 80*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 2)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) - 5)*(x^2/4 - x/4 + y^2/4 - y/4 + 1/8))/100000000 + ((x^2 - x + y^2 - y + 1/2)^2*(4*(3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 1)*(2*x^2 - 2*x + 2*y^2 - 2*y + 1) - (3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 + 80*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 2)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) - 1)*(x^2/2 - x/2 + y^2/2 - y/2 + 1/4) - (tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2))/20 + 160*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 2)*(2*x^2 - 2*x + 2*y^2 - 2*y + 1)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 1/800))/2 + ((3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 1)*(2*x^2 - 2*x + 2*y^2 - 2*y + 1))/12500000 + (tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2))/500000000 - (3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2)/2000000000000 + (pi*sin((pi*t)/3)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2)*(x^2 - x + y^2 - y + 5001/10000)^3)/960 - 99/2000000000000))/(3*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2)*(2*x^2 - 2*x + 2*y^2 - 2*y + 5001/5000)^3)",2);
    // gsFunctionExpr<> sourceQnew("(40*(192*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12) - 192)*(tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12) + 1)*((3*(4*(3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 1)*(5*x^2 - 5*x + 5*y^2 - 5*y + 5/2) - (3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 + 40*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 2)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) - 1)*(x^2 - x + y^2 - y + 1/2) - (tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2))/5 + 80*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 2)*(2*x^2 - 2*x + 2*y^2 - 2*y + 1)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 1/200)*(x^2/3 - x/3 + y^2/3 - y/3 + 1/6))/10000 - ((15*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 + 80*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 2)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) - 5)*(x^2/4 - x/4 + y^2/4 - y/4 + 1/8))/100000000 + ((x^2 - x + y^2 - y + 1/2)^2*(4*(3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 1)*(2*x^2 - 2*x + 2*y^2 - 2*y + 1) - (3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 + 80*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 2)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) - 1)*(x^2/2 - x/2 + y^2/2 - y/2 + 1/4) - (tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2))/20 + 160*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 2)*(2*x^2 - 2*x + 2*y^2 - 2*y + 1)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 1/800))/2 + ((3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 1)*(2*x^2 - 2*x + 2*y^2 - 2*y + 1))/12500000 + (tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2))/500000000 - (3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2)/2000000000000 + (pi*sin((pi*t)/3)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2)*(x^2 - x + y^2 - y + 5001/10000)^3)/960 - 99/2000000000000))/(3*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2)*(2*x^2 - 2*x + 2*y^2 - 2*y + 5001/5000)^3)",2);
    gsFunctionExpr<> sourceQold("(40*(192*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12) - 192)*(tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12) + 1)*((3*(((3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 1)*(5*x^2 - 5*x + 5*y^2 - 5*y + 5/2))/2 - (3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 + 40*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 2)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) - 1)*(x^2 - x + y^2 - y + 1/2) - (tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2))/40 + 10*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 2)*(2*x^2 - 2*x + 2*y^2 - 2*y + 1)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 1/1600)*(x^2/3 - x/3 + y^2/3 - y/3 + 1/6))/10000 - ((15*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 + 80*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 2)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) - 5)*(x^2/4 - x/4 + y^2/4 - y/4 + 1/8))/100000000 + ((x^2 - x + y^2 - y + 1/2)^2*(((3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 1)*(2*x^2 - 2*x + 2*y^2 - 2*y + 1))/2 - (3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 + 80*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 2)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) - 1)*(x^2/2 - x/2 + y^2/2 - y/2 + 1/4) - (tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2))/160 + 20*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 2)*(2*x^2 - 2*x + 2*y^2 - 2*y + 1)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 1/6400))/2 + ((3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 1)*(2*x^2 - 2*x + 2*y^2 - 2*y + 1))/100000000 + (tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2))/4000000000 - (3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2)/2000000000000 + (pi*sin((pi*t)/3)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2)*(x^2 - x + y^2 - y + 5001/10000)^3)/960 - 23/4000000000000))/(3*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2)*(2*x^2 - 2*x + 2*y^2 - 2*y + 5001/5000)^3)",2);
    gsFunctionExpr<> sourceQnew("(40*(192*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12) - 192)*(tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12) + 1)*((3*(((3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 1)*(5*x^2 - 5*x + 5*y^2 - 5*y + 5/2))/2 - (3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 + 40*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 2)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) - 1)*(x^2 - x + y^2 - y + 1/2) - (tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2))/40 + 10*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 2)*(2*x^2 - 2*x + 2*y^2 - 2*y + 1)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 1/1600)*(x^2/3 - x/3 + y^2/3 - y/3 + 1/6))/10000 - ((15*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 + 80*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 2)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) - 5)*(x^2/4 - x/4 + y^2/4 - y/4 + 1/8))/100000000 + ((x^2 - x + y^2 - y + 1/2)^2*(((3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 1)*(2*x^2 - 2*x + 2*y^2 - 2*y + 1))/2 - (3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 + 80*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 2)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) - 1)*(x^2/2 - x/2 + y^2/2 - y/2 + 1/4) - (tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2))/160 + 20*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 2)*(2*x^2 - 2*x + 2*y^2 - 2*y + 1)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 1/6400))/2 + ((3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2 - 1)*(2*x^2 - 2*x + 2*y^2 - 2*y + 1))/100000000 + (tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2))/4000000000 - (3*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2)/2000000000000 + (pi*sin((pi*t)/3)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2)*(x^2 - x + y^2 - y + 5001/10000)^3)/960 - 23/4000000000000))/(3*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2)*(2*x^2 - 2*x + 2*y^2 - 2*y + 5001/5000)^3)",2);

    gsVector<> pt2(2,1); pt2<<0.5, 0.5;
    sourceQold.set_t(0);
    auto Q_eval= ev.getVariable(sourceQold, G);
    gsInfo<< "Value at of Q :" << ev.eval(Q_eval, pt2) <<"\n";

    // Derivatives of the polynomial double well potential (M. Kästner et al., 2016)
    auto dmu_c = - 1.0 + 3.0 * (c*c).val(); // f_2 (second derivative of double well)
    auto ddmu_c = 6*c.val(); // f_3 (third derivative of double well)

    // Mobility
    auto M_c  = 1.0 + 0.0*c.val(); // replace with const_expr(1.0) instead of using 0*c
    auto dM_c = 0.0 * gradc; // replace with const_expr(1.0) instead of using 0*c!!

    auto residual = w*dc + // M
                M_c.val() * igrad(w,G) * dmu_c * gradc.tr() + // F_bar
                M_c.val() * ilapl(w,G) *lambda *laplc.val(); // K_laplacian

    //  =========== Terms for boundary integrals ===========
    // (1) Neumann boundary condition
    gsFunctionExpr<> bc1("-(320*(x - 1/2)*(tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12) - 1)*(tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12) + 1)*((80033*y)/80000 + (x^2 - x + y^2 - y + 5001/10000)^2 + (x*(6400*y^2 - 6400*y + 80033/25))/3200 - (x^2*(6400*y^2 - 6400*y + 160033/25))/3200 - (tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2*(6*x^2 - 6*x + 6*y^2 - 6*y + 15003/5000))/20000 + 2*x^3 - x^4 - (160033*y^2)/80000 + 2*y^3 - y^4 + (tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2)*(x^2 - x + y^2 - y + 1251/2500))/80 - 2001651/8000000))/((4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2)*(2*x^2 - 2*x + 2*y^2 - 2*y + 5001/5000)^2)","-(320*(y - 1/2)*(tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12) - 1)*(tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12) + 1)*((80033*x)/80000 + (x^2 - x + y^2 - y + 5001/10000)^2 + (y*(6400*x^2 - 6400*x + 80033/25))/3200 - (y^2*(6400*x^2 - 6400*x + 160033/25))/3200 - (tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)^2*(6*x^2 - 6*x + 6*y^2 - 6*y + 15003/5000))/20000 - (160033*x^2)/80000 + 2*x^3 - x^4 + 2*y^3 - y^4 + (tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2)*(x^2 - x + y^2 - y + 1251/2500))/80 - 2001651/8000000))/((4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2)*(2*x^2 - 2*x + 2*y^2 - 2*y + 5001/5000)^2)",2);

    // (2) Laplace boundary condition
    gsFunctionExpr<> bc2("(160*(tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12) - 1)*(tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12) + 1)*(x^2 - y - x + y^2 + 40*tanh(cos((pi*t)/3) - 20*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2) + 12)*(4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2)*(x^2 - x + y^2 - y + 1/2) + 2501/5000))/((4*x^2 - 4*x + 4*y^2 - 4*y + 5001/2500)^(1/2)*(2*x^2 - 2*x + 2*y^2 - 2*y + 5001/5000))",2);
    // =====================================================

    // ![Initialize the assembler]
    w.setup(bc, dirichlet::l2Projection, 0);
    A.initSystem();
    // ![Initialize the assembler]

    // Define linear solver (install SuperLUMT-devel)
#ifdef GISMO_WITH_SUPERLU
    gsSparseSolver<>::SuperLU solver;
#   else
    gsSparseSolver<>::LU solver;
#endif

    gsMatrix<> Q;
    gsSparseMatrix<> K;

    // Legend:
    // C_old   = C_n
    // C       = C_n+1,i-1
    // C_alpha = C_{n+alpha_f,i-1}
    // dC

    gsInfo<<"Starting.."<<"\n";
    GISMO_ASSERT(projection_Crs<=2,"Projection method not implemented, index should be 0 (L2) or 1 (Local QI) for coarsening, but is "<<projection_Crs);

    gsInfo<<"Initial condition.."<<"\n";

    if (random)
    {
        // // %%%%%%%%%%%%%%%%%%%%%%%% Random initial condition %%%%%%%%%%%%%%%%%%%%%%%%
        // gsMatrix<> tmp = gsMatrix<>::Random(A.numDofs(),1);
        // // gsMatrix<> tmp = gsMatrix<>::Random(dbasis.size(),1);
        // tmp.array() *= CHopt.askReal("ampl",0.005); //random uniform variable in [-0.05,0.05]
        // tmp.array() += CHopt.askReal("mean",0.0); // 0.45
        // mp_cold.addPatch(dbasis.basis(0).makeGeometry(tmp));
        // tmp.setZero();
        // mp_dcold.addPatch(dbasis.basis(0).makeGeometry(tmp));
        // gsInfo<<mp_cold.patch(0).coefs().size()<<"\n";
        // gsInfo<<mp_dcold.patch(0).coefs().size()<<"\n";

        // // %%%%%%%%%%%%%%%%%%%%%%%% XML initial condition %%%%%%%%%%%%%%%%%%%%%%%%
        // gsFileData<> fd1;
        // std::string file_name;
        // if (pattern==0) // nucleation
        //     file_name = "/Users/lucasventavinuela/gismo/build/nucleation_cl1024.xml";
        // else
        //     file_name = "/Users/lucasventavinuela/gismo/build/new_spin.xml";

        // gsMultiBasis<> dbasis_IC;
        // gsMatrix<> Coefs;

        // fd1.read(file_name);
        // fd1.getId(0,dbasis_IC);
        // fd1.getId(1,Coefs);

        // gsGeometry<>::uPtr IC_function = dbasis_IC.basis(0).makeGeometry(give(Coefs));
        // gsQuasiInterpolate<real_t>::localIntpl(dbasis.basis(0),*IC_function,Cold);

        // // // New: to update the values of the free DoFs
        // // const gsDofMapper & mapper = w.mapper();
        // // gsGeometry<>::uPtr geom  = dbasis.basis(0).makeGeometry(give(Cold));
        // // Cold.resize(mapper.freeSize(),1);

        // // for (index_t c = 0; c!=geom->coefs().cols(); c++) // for all components
        // // {
        // //     for (index_t i = 0; i != geom->coefs().rows(); ++i)
        // //     {
        // //         const index_t ii = mapper.index(i, c);
        // //         if ( mapper.is_free_index(ii) ) // DoF value is in the solVector
        // //             Cold.at(ii) = geom->coefs()(i, c);
        // //     }
        // // }

        // gsInfo<<"COLD0"<<Cold<<"\n";
        gsFileData<> fd1;
        std::string file_name;
        if (pattern==0) // nucleation
            file_name = "multipatch.xml";
        else
            file_name = "/Users/lucasventavinuela/gismo/build/new_spin.xml";

       gsMultiPatch<> MP;

        fd1.read(file_name);
        fd1.getId(0,MP);
        gsInfo<<"Imported multipatch\n";

        mp_cold.addPatch(MP.patch(0));
        Cold.setZero(dbasis.size(),1);
        mp_dcold.addPatch(dbasis.basis(0).makeGeometry(Cold));

        // gsWriteParaview(mp,mp_cold, out+"/initial_condition", 5000);

    }
    else
    {
        // %%%%%%%%%%%%%%%%%%%%%%%% Analytical intial condition %%%%%%%%%%%%%%%%%%%%%%%%
        GISMO_ASSERT(mp.geoDim()==source.domainDim(),"Domain dimension of the source function should be equal to the geometry dimension, but "<<source.domainDim()<<"!="<<mp.geoDim());
        gsMatrix<> tmp;
        Cold.setZero(A.numDofs(),1);
        gsFunctionExpr<> initial_cond("tanh(cos((pi*t)/3) - 40*((x - 1/2)^2 + (y - 1/2)^2 + 1/10000)^(1/2) + 12)",2);
        initial_cond.set_t(0);
        real_t error = gsL2Projection<real_t>::project(dbasis,mp,initial_cond,tmp,A.options());  // 3rd arg has to be multipatch
        if (verbose>0) gsInfo << "L2 projections error "<<error<<"\n";
        mp_cold.addPatch(dbasis.basis(0).makeGeometry(tmp));
        tmp.setZero();
        mp_dcold.addPatch(dbasis.basis(0).makeGeometry(tmp));
    }

    gsWriteParaview(mp,mp_cold,out+"/initial_condition",5000);

    mp_cnew = mp_cold;
    mp_dcnew = mp_dcold;

    // Fill multipatch zeros!
    mp_qnew = mp_dcold;
    mp_qold = mp_dcold;

    real_t Q0norm = 1, Qnorm = 10;
    real_t tol = TIMEopt.askReal("tol",1e-4);

    gsParaviewCollection collection(out+"/solution", &ev);
    collection.options().setSwitch("plotElements", true);
    collection.options().setInt("plotElements.resolution", 4);
    collection.options().setInt("numPoints",(mp.geoDim()==3) ? 10000 : 5000);
    collection.options().setInt("precision", 10); // digits 10^-10

    gsParaviewCollection collection_mesh(out+"/collect", &ev);
    collection_mesh.options().setSwitch("plotElements", true);
    collection_mesh.options().setInt("plotElements.resolution", 4);
    collection_mesh.options().setInt("numPoints",(mp.geoDim()==3) ? 10000 : 5000);
    collection_mesh.options().setInt("precision", 10); // digits 10^-10
    // new collection for errors
    // gsParaviewCollection error_collection("ParaviewOutput/errors", &ev);
    // error_collection.options().setSwitch("plotElements", true);
    // error_collection.options().setInt("plotElements.resolution", 4);
    // error_collection.options().setInt("numPoints",(mp.geoDim()==3) ? 10000 : 5000);

    real_t time = 0;
    bool converged = false;

    // Sparse matrix for Nitsche contribution
    gsSparseMatrix<> K_nitsche; // empty variable


    // Errors projection
    real_t error_crs_c = 0;
    real_t error_crs_dc = 0;

    real_t error_ref_cnew = 0;
    real_t error_ref_dcnew = 0;
    real_t error_ref_cold = 0;
    real_t error_ref_dcold = 0;

    // ! [Load mesher options]
    // DIMENSION-INDEPENDENT
    // gsAdaptiveMeshingBase<real_t> * mesher;
    // if(mp.geoDim()==2)
    //     mesher = new gsAdaptiveMeshing<dim,real_t>(dbasis);
    // else if(mp.geoDim()==2)
    //     mesher = new gsAdaptiveMeshing<dim,real_t>(dbasis);

    gsAdaptiveMeshing<dim,real_t> mesher;
    if (MESHopt.askSwitch("Adaptive",true))
    {
        mesher = gsAdaptiveMeshing<dim,real_t>(dbasis);
        mesher.options().setInt("RefineRule",MESHopt.askInt("RefineRule",1));
        mesher.options().setInt("CoarsenRule",MESHopt.askInt("CoarsenRule",1));
        mesher.options().setReal("RefineParam",MESHopt.askReal("RefineParam",0.1));
        mesher.options().setReal("CoarsenParam",MESHopt.askReal("CoarsenParam",0.1));
        mesher.options().setSwitch("Admissible",MESHopt.askSwitch("Admissible",false));
        mesher.options().setInt("Jump",MESHopt.askInt("Jump",2));
        mesher.options().setInt("MaxLevel",numRefine);
        mesher.options().setSwitch("Absolute",MESHopt.askSwitch("Absolute",true));
        mesher.getOptions();
    }
    // ! [Load mesher options]

    std::ofstream csvFile;
    csvFile.open(out+"/out.csv");
    csvFile << std::left
            << std::setw(12) << "TimeStep"        << " , "
            << std::setw(10) << "RefIt"           << " , "
            << std::setw(12) << "NumDOFs"         << " , "
            << std::setw(14) << "Mass"            << " , "
            << std::setw(16) << "nSolvesStep"     << " , "
            << std::setw(20) << "l2err - conv"    << " , "
            << std::setw(20) << "h1err - conv"    << " , "
            << std::setw(20) << "h2err - conv"    << " , "
            << std::setw(15) << "proj_err c"      << " , "
            << std::setw(15) << "proj_err dc"     << " , "
            << std::setw(12) << "AT"              << " , "
            << std::setw(12) << "ST"              << " , "
            << "PT\n";

    gsVector<> pt(2,1); pt<<0.5, 0.5;

    // auto u_manufactured = ev.getVariable(source_time, G);
    // gsInfo<< "Value at mid-point :" << ev.eval(neumann_manuf, pt) <<"\n";

    real_t mass2 = ev.integral(meas(G)*cold);
    gsInfo<<mass2<<"\n";
    real_t told, tnew;

    real_t assemblyTimestep, solverTimestep, projTimestep;
    real_t assemblyTimeRefIt, solverTimeRefIt, projTimeRefIt;

    gsFunctionExpr<> source_time("tanh(cos((pi*t)/3) - 40*((x - 1/2)^2 + (y - 1/2)^2 + 1/10000)^(1/2) + 12)",2);

    for (index_t step = 0; step!=maxSteps; step++)
    {   

        gsDebugVar(mp_cold.patch(0).coefs().maxCoeff());
        gsDebugVar(mp_cold.patch(0).coefs().minCoeff());
        
        assemblyTimestep = 0;
        solverTimestep   = 0;
        projTimestep     = 0;
        
        for (index_t refIt = 0; refIt!=MESHopt.askInt("RefIt",5); refIt++)
        {
            nSolvesStep = 0;
            assemblyTimeRefIt = 0;
            solverTimeRefIt   = 0;
            projTimeRefIt = 0;

            if (MESHopt.askSwitch("Adaptive",true))
            {
                clock.restart();
                if (projection_Crs == 0)
                {
                    error_crs_c = 0;
                    error_crs_dc = 0;
                    error_crs_c =  gsL2Projection<real_t>::project(dbasis, ibasis, mp, mp_cold.patch(0), CnewF, A.options());
                    error_crs_dc = gsL2Projection<real_t>::project(dbasis, ibasis, mp, mp_dcold.patch(0), dCnewF, A.options());   
                    gsDebug<<"Error in the L2 projection of the initial condition (c) : "<<error_crs_c<<"\n";
                    gsDebug<<"Error in the L2 projection of the initial condition (dc): "<<error_crs_dc<<"\n";  
                    // gsWriteParaview(mp,mp_cold, out+"/step_"+std::to_string(step)+"_refit_"+std::to_string(refIt), 5000);
                    // gsWriteParaview(mp_cold.,dbasis,5000);
                }
                else if (projection_Crs == 1)
                {
                    gsQuasiInterpolate<real_t>::localIntpl(dbasis.basis(0),mp_cold.patch(0),CnewF);
                    gsQuasiInterpolate<real_t>::localIntpl(dbasis.basis(0),mp_dcold.patch(0),dCnewF);    
                }
                projTimeRefIt += clock.stop();
            }
            else
            {
                // Tensor product solution!
                CnewF = mp_cold.patch(0).coefs();
                dCnewF = mp_dcold.patch(0).coefs();
            }

            // sourceQ.set_t(step*dt); //update t in the source term
            // auto     qq = A.getCoeff(sourceQ,G);
                      
            // Resize the data structure inside the mesher
            if (MESHopt.askSwitch("Adaptive",true))
                mesher.rebuild();

            // Setup the space (compute Dirichlet BCs)
            w.setup(bc, dirichlet::l2Projection, 0);

            // Reset the assembler
            A.initSystem();
            //A.initVector();

     

            for (index_t dt_it = 0; dt_it != lmax; dt_it++)
            {
                if (verbose>0) gsInfo<<"Time step "<<step<<"/"<<maxSteps<<", iteration "<<dt_it<<": dt = "<<dt<<", [t_start,t_end] = ["<<time<<" , "<<time+dt<<"]"<<"\n";
                tmp_alpha_m = tmp_alpha_f = tmp_gamma = 1;
                tmp_alpha_m_func.setValue(tmp_alpha_m,dim);
                tmp_alpha_f_func.setValue(tmp_alpha_f,dim);


                for (index_t k = 0; k!=2; k++)
                {
                    converged = false;
                    std::string method = (k==0) ? "Backward Euler " : "Generalized Alpha ";

                    Cnew.setZero(A.numDofs(),1);
                    dCnew.setZero(A.numDofs(),1);

                    for (index_t i = 0; i < dbasis.basis(0).size(); i++)
                    {
                        if (w.mapper().is_free(i))
                        {
                            Cnew(w.mapper().index(i),0) = CnewF(i,0);
                            dCnew(w.mapper().index(i),0) = (tmp_gamma-1)/tmp_gamma * dCnewF(i,0);
                        }
                    }

                    //  ////// ===========
                    //  source_time.set_t(tnew);
                    //  auto u_manufactured = ev.getVariable(source_time, G);
                    //  real_t l2err3;
                    //  l2err3 = math::sqrt( ev.integral( (u_manufactured - cnew_sol).sqNorm() * meas(G) ) ); // / ev.integral(ff.sqNorm()*meas(G)) );
                    //  gsInfo << "L2 error after refinement: " << l2err3 << "\n";
                    //  // ==========================
 
                    
                    Q0norm = 1;
                    Qnorm = 10;

                    // ==== For Nitsche BC ====
                    Cold.setZero(A.numDofs(),1);
                    dCold.setZero(A.numDofs(),1);
                    Cold = Cnew;
                    dCold = dCnew;
                    // ========================

                    // ====== Time integration of the source term ======
                    told = time;
                    tnew = time + dt;
                    real_t alpha_time = (1-tmp_alpha_f)*told + tmp_alpha_f*tnew;
                    sourceQold.set_t(told);
                    sourceQnew.set_t(tnew);
                    auto     qalpha = (1-tmp_alpha_f)*A.getCoeff(sourceQold,G) + tmp_alpha_f*A.getCoeff(sourceQnew,G);
                    // =================================================

                    for (index_t it = 0; it!= maxIt; it++)
                    {
                        A.initMatrix();
                        A.initVector();
                        // A.clearRhs(); // Resets to zero the values of the already allocated to residual (RHS)
                        A.clearRhs(); // Resets to zero the values of the already allocated to residual (RHS)
                        // ==== For Nitsche BC ====
                        Calpha.noalias()  = Cold  + tmp_alpha_f * ( Cnew  - Cold );
                        dCalpha.noalias() = dCold + tmp_alpha_m * ( dCnew - dCold);
                        // ========================

                        clock.restart();
                        // Assemble the RHS
                        A.assemble(residual * meas(G));
                        
                        // ==== Assemble source term (manufactured) ====
                        A.assemble(meas(G)* (- w) * qalpha);
                        // =============================================
                        assemblyTimeRefIt += clock.stop();

                        Q = A.rhs();

                        // ==== Assemble boundary terms (manufactured) ====
                        // Evaluate the boundary terms at alpha time 
                        clock.restart();
                        bc1.set_t(alpha_time);
                        auto bc_Neumann = A.getCoeff(bc1,G); // Neumann BC
                        A.assembleBdr(bc.get("Neumann"), w * (bc_Neumann.tr() * nv(G)));  // assemble boundary term -- flux from manufactured solutions
                        
                        bc2.set_t(alpha_time);
                        auto bc_Laplace = A.getCoeff(bc2,G); // Neumann BC
                        A.assembleBdr(bc.get("Laplace"), (igrad(w,G).tr() * nv(G))* lambda * bc_Laplace.tr());
                        assemblyTimeRefIt += clock.stop();
                        // ================================================

                        // // Assemble the Nitsche BC on the sides with Neumann condition
                        // // A.clearMatrix(); // Resets to zero the values of the already allocated to matrix (LHS)
                        // A.initMatrix();
                        // clock.restart();
                        // A.assembleBdr(bc.get("Neumann"), - lambda * igrad(w,G) *  nv(G)  * ilapl(w,G).tr() + // consistency term
                        //             penalty * (igrad(w,G) * nv(G).normalized()) * hmax * (igrad(w,G) * nv(G)).tr() - // penalty (stabilizing) term
                        //             lambda * ilapl(w,G) * (igrad(w,G)  * nv(G)).tr()); // symmetry term            

                        // assemblyTime += clock.stop();
                        // K_nitsche = A.giveMatrix(); // .giveMatrix() moves the matrix A into K_nitche (avoids having two matrices A and K_nitsche)

                        // if (bc.get("Neumann").size()!=0)
                        //     Q.noalias() += K_nitsche * Calpha; // add the residual term from Nitche (using the matrix )

                        A.initMatrix();
                        // Assembly of the tangent stiffness matrix (K_m and K_f simultaneously) %%
                        clock.restart();
                        A.assemble(meas(G) * (w*w.tr()*tmp_alpha_m +// K_m
                                            (tmp_alpha_f * tmp_gamma * dt)* (dmu_c *igrad(w,G) * igrad(w,G).tr() + // K_f1
                                            ddmu_c * igrad(w,G) * gradc.tr() * w.tr() + // K_f2
                                            lambda * ilapl(w,G) * ilapl(w,G).tr()))); // K_laplacian
                                            // lambda * igrad(w,G)*dM_c.tr()*ilapl(w,G).tr()   +  // K_mobility
                        assemblyTimeRefIt += clock.stop();

                        K = A.giveMatrix();
                        // if (bc.get("Neumann").size()!=0)
                        //     K += (tmp_alpha_f * tmp_gamma * dt) * K_nitsche; // add the Nitsche term to the stiffness matrix


                        clock.restart();
    #ifdef GISMO_WITH_SUPERLU
                        if (0==k)
                            solver.analyzePattern(K);
                        solver.factorize(K);
    #else
                        solver.compute(K);
    #endif

    // #ifdef GISMO_WITH_PARDISO
    // gsSparseSolver<>::PardisoLU solverpLU;

                        dCupdate = solver.solve(-Q);
                        solverTimeRefIt += clock.stop();
                        nSolves++;
                        nSolvesStep++;

                        dCnew += dCupdate;
                        Cnew.noalias() += (tmp_gamma*dt)*dCupdate;

                        // ====================== CHECK THE CONVERGENCE CONDITIONS ======================
                        if (it == 0) Q0norm = Q.norm();
                        else         Qnorm = Q.norm();

                        real_t dCnorm = dCupdate.norm();
                        real_t Cnewnorm = Cnew.norm();
                    
                        if (verbose==2) gsInfo<<"\t\tNR iter   "<<it<<": res  = "<<Qnorm/Q0norm<<"\n";
                        if (verbose==2) gsInfo<<"\t\t          "<<it<<": dc/c = "<<dCnorm/Cnewnorm<<"\n";

                        if (it>0 && Qnorm/Q0norm < tol && dCnorm/Cnewnorm < tol)
                        {
                            if (verbose>0) gsInfo<<"\t\t"<<method<<"converged in "<<it<<" iterations\n";
                                converged = true;
                            break;
                        }
                        else if (it==maxIt-1)
                        {
                            if (verbose>0) gsInfo<<"\t\t"<<method<<"did not converge!\n";
                                converged = false;
                            break;
                        }
                        // ===============================================================
                    }
                    if (!converged)
                        break;

                    // %% Switch to generalized-alpha parameters (k=1)
                    tmp_alpha_m = alpha_m;
                    tmp_alpha_f = alpha_f;
                    tmp_alpha_m_func.setValue(tmp_alpha_m,dim);
                    tmp_alpha_f_func.setValue(tmp_alpha_f,dim);
                    tmp_gamma = gamma;

                    // %% For time step adaptivity %%
                    // Csols[k] = Cnew; // k=0: BE, k=1: alpha
                }// Backward Euler/Generalized Alpha

                // %%%%%%%%%% For time step adaptivity %%%%%%%%%%
                // if (converged)
                // {
                //     t_err = (Csols[0] - Csols[1]).norm() / (Csols[1]).norm();
                //     dt_old = dt;
                //     // dt *= t_rho * math::sqrt(TOL / t_err);
                //     if (t_err < TOL)
                //         break;
                // }
                // else
                // {
                //     dt_old = dt;
                //     // dt *= t_rho;
                // }
            }// time step adaptivity

            // Update the old c and dc into splines
            cnew_sol.extract(mp_cnew);
            dcnew_sol.extract(mp_dcnew);

            // mp_qnew = mp_qold;

            real_t mass = ev.integral(meas(G)*cnew);
            source_time.set_t(tnew);
            auto u_manufactured = ev.getVariable(source_time, G);
           //  gsWriteParaview(source_time,out+"/manu_sol");
            // ==========================
            real_t l2err2, h1err2, h2err2;
            l2err2 = math::sqrt( ev.integral( (u_manufactured - cnew_sol).sqNorm() * meas(G) ) ); // / ev.integral(ff.sqNorm()*meas(G)) );
            h1err2 = l2err2 + math::sqrt(ev.integral((igrad(u_manufactured) - igrad(cnew_sol, G)).sqNorm() * meas(G)));
            h2err2 = h1err2 + math::sqrt(ev.integral( ( ihess(u_manufactured) - ihess(cnew_sol,G) ).sqNorm() * meas(G) )); // /ev.integral( ihess(ff).sqNorm()*meas(G) )

            csvFile << std::left
                    << std::setw(12) << step              << " , "
                    << std::setw(10) << refIt             << " , "
                    << std::setw(12) << A.numDofs()       << " , "
                    << std::setw(14) << mass              << " , "
                    << std::setw(16) << nSolvesStep       << " , "
                    << std::setw(20) << l2err2            << " , "
                    << std::setw(20) << h1err2            << " , "
                    << std::setw(20) << h2err2            << " , "
                    << std::setw(15) << error_crs_c       << " , "
                    << std::setw(15) << error_crs_dc      << " , "
                    << std::setw(12) << assemblyTimeRefIt  << " , "
                    << std::setw(12) << solverTimeRefIt    << " , "
                    << projTimeRefIt << "\n";

            csvFile.flush();  // Forces the file to write immediately
            
            solverTimestep     += solverTimeRefIt;
            assemblyTimestep   += assemblyTimeRefIt;
            projTimestep       += projTimeRefIt;
            gsInfo << "Time step " << step << " took pro " <<  projTimeRefIt << " seconds.\n";
            gsInfo << "Time step " << step << " took ass " <<  assemblyTimeRefIt << " seconds.\n";
            gsInfo << "Time step " << step << " took sol " <<  solverTimeRefIt << " seconds.\n"; 

            // Plot the mesh and solution at each refinement iteration
            // I am plotting the new solution on the mesh for a refIt!
            if (plot)
            {
                // Export the mesh
                collection_mesh.newTimeStep(&mp);
                auto multipatch_projected = ev.getVariable(mp_cnew, G);
                collection_mesh.addField(multipatch_projected,"solu");
                collection_mesh.saveTimeStep();
            }


            if (MESHopt.askSwitch("Adaptive",true))
            {
                // -------------REFINEMENT-------------------
                // Compute the integral of c over each element
                ev.integralElWise(meas(G) * cnew_sol);
                std::vector<real_t> cInt = ev.elementwise();
                gsAsVector<real_t> cvec(cInt.data(),cInt.size());  // Temporary Eigen::Map
                // Compute the area of each element
                ev.integralElWise(meas(G));
                std::vector<real_t> areas = ev.elementwise();
                gsAsVector<real_t> avec(areas.data(),areas.size()); // Temporary Eigen::Map

                // Invert and normalize the element-wise average (c/area), as:
                // err = 1-|c|/a;
                cvec.array() = 1-(cvec.array().abs()/avec.array());

                // Mark the elements for refinement
                gsHBoxContainer<dim,real_t> refine;//, coarsen; // llamar markedRef???????? en vez de refine
                mesher.markRef_into(cInt,refine);

                // If elements are marked for refinement
                if (refine.totalSize()!=0)
                {
                    // Refine dbasis
                    if (verbose>1) gsInfo<<"Basis before refinement:\n "<<dbasis.basis(0)<<"\n";
                    mesher.refine(refine);
                    if (verbose>1) gsInfo<<"Basis after refinement:\n "<<dbasis.basis(0)<<"\n";
                }// refine
                else
                    break;
            } // refinement switch
            else
                break; // break the refinement loop
        }// mesh adaptivity

        //! [Export visualization in ParaView]
        if (plot && step % plotmod==0)
        {
            // Export the mesh
            collection.newTimeStep(&mp);
            collection.addField(cnew,"numerical solution");
            source_time.set_t(tnew);
            auto u_manufactured = ev.getVariable(source_time, G);
            collection.addField(u_manufactured,"analytical solution");
            auto mp_var = ev.getVariable(mp_cnew, G);
            collection.addField((u_manufactured - mp_var).sqNorm(), "L2 error");
            collection.addField((igrad(u_manufactured) - igrad(mp_var, G)).sqNorm(), "H1 seminorm error");
            collection.addField((u_manufactured - mp_var).sqNorm() + (igrad(u_manufactured) - igrad(mp_var, G)).sqNorm(), "H1 full error");
            gsInfo << "Number of degrees of freedom:\t" << A.numDofs()  << std::endl;
            collection.saveTimeStep();
        }

        // real_t mass = ev.integral(meas(G)*cnew);
        // csvFile << step << "," << A.numDofs() <<"," << mass <<  "," << error_ref_cnew << ","<< error_ref_dcnew << ","<<error_crs_c<<","<<error_crs_dc<< "\n";
        // // Reset the errors after the loop
        // error_ref_cnew  = 0;
        // error_ref_cold  = 0;
        // error_ref_dcnew = 0;
        // error_ref_dcold = 0;

        // mp_cold = mp_cnew; // mp_cold has the basis of mp_cnew (fine basis)
        // mp_dcold = mp_dcnew; // mp_dcold has the basis of mp_dcnew (fine basis)

        if (MESHopt.askSwitch("Adaptive",true) && step!=maxSteps)
        {
            // The mesh changed due to refinement, so w changed.
            // Hence we need to call setup
            // w.setup(bc, dirichlet::l2Projection, 0);
            // // In addition, since the size of the space changes,
            // // the size of the corresponding solution vector
            // // is wrong. Therefore, we need a projection,
            // // which is exact since it is on a finer mesh.
            // // We only need to do this for the solution vector,
            // // since that is the one that we use to compute the errors
            // // for coarsening.
            // real_t ref_proj = interpolate_ref(dbasis, mp_cnew.patch(0), CnewF, mp);
            // // w.setup(bc, dirichlet::l2Projection, 0);

            // Cnew.resize(w.mapper().freeSize(),1); // freeSize if updated because of the setup
            // // Insert the interpolated coefficients inside the solution objects
            // for (index_t i = 0; i < dbasis.basis(0).size(); i++)
            //     if (w.mapper().is_free(i))
            //         Cnew(w.mapper().index(i),0) = CnewF(i,0);


            // If the mesh adaptivity is converged, we perform a coarsening step
            // -------------COARSENING-------------------
            // Resize the mesher data structure
            mesher.rebuild();

            // Compute the integral of c over each element
            ev.integralElWise(meas(G) * cnew);
            std::vector<real_t> cInt = ev.elementwise();
            gsAsVector<real_t> cvec(cInt.data(),cInt.size());  // Temporary Eigen::Map
            // Compute the area of each element
            ev.integralElWise(meas(G));
            std::vector<real_t> areas = ev.elementwise();
            gsAsVector<real_t> avec(areas.data(),areas.size()); // Temporary Eigen::Map

            // Invert and normalize the element-wise average (c/area), as:
            // err = 1-|c|/a;
            cvec.array() = 1-(cvec.array().abs()/avec.array());

            // Coarsen everything above threshold (opposite of refinement)
            gsHBoxContainer<dim,real_t> coarsen;
            mesher.markCrs_into(cInt,coarsen); // includes admissibility

            // If elements are marked for refinement
            if (coarsen.totalSize()!=0)
            {
                // Refine dbasis
                if (verbose>1) gsInfo<<"Basis before coarsening:\n "<<dbasis.basis(0)<<"\n";
                mesher.unrefine(coarsen);
                if (verbose>1) gsInfo<<"Basis after coarsening:\n "<<dbasis.basis(0)<<"\n";
            }// coarsen
        } // coarsening switch

        // Update time and old solutions
        solverTime     += solverTimestep;
        assemblyTime   += assemblyTimestep;
        projectionTime += projTimestep;


        time += dt_old;
        mp_cold = mp_cnew;
        mp_dcold = mp_dcnew;

        // gsDebugVar(mp_cold.patch(0).coefs().maxCoeff());
        // gsDebugVar(mp_cold.patch(0).coefs().minCoeff());
    }
    if (plot)
    {
        collection.save();
        collection_mesh.save();
    }
    // else if (plot_error)
    //     error_collection.save();
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";


    csvFile << "TOTAL COMPUTATIONAL TIMES, Value\n";
    csvFile << "Assembly  :" << assemblyTime << "\n";
    csvFile << "Solver    :" << solverTime << "\n";
    csvFile << "Projection:" << projectionTime << "\n";
    csvFile << "Total     :" << assemblyTime + solverTime + projectionTime << "\n";
    csvFile << "NSolver   :" << nSolves << "\n";

    gsInfo<<"[CLOCK] --- Time for assembly: "<<assemblyTime<<" [s]\n";
    gsInfo<<"[CLOCK] --- Time for solver  : "<<solverTime<<" [s]\n";
    gsInfo<<"[CLOCK] --- Time for projection: "<<projectionTime<<" [s]\n";
    gsInfo<<"[CLOCK] --- Number of solves: "<<nSolves<<"\n";
  
     real_t mass = ev.integral(meas(G)*cnew);
 
 
     // ========= Compute L2 error ========= 
     // === Solution with time ===
     source_time.set_t(time);
     gsDebugVar(time);
     auto u_manufactured = ev.getVariable(source_time, G);
     auto mp_var = ev.getVariable(mp_cnew, G);
     // ==========================
     real_t l2err, h1err, h2err;     
     l2err  = math::sqrt( ev.integral( (u_manufactured - mp_var).sqNorm() * meas(G) ) ); // / ev.integral(ff.sqNorm()*meas(G)) );
     h1err = l2err + math::sqrt(ev.integral((igrad(u_manufactured) - igrad(mp_var, G)).sqNorm() * meas(G)));
     h2err = h1err + math::sqrt(ev.integral((ihess(u_manufactured) - ihess(mp_var, G)).sqNorm() * meas(G) )); // /ev.integral( ihess(ff).sqNorm()*meas(G) )
     gsInfo << " L2 error : " << l2err << "\n";
     gsInfo << " H1 error : " << h1err << "\n";
     gsInfo << " H2 error : " << h2err << "\n";
     gsInfo << " Mesh-size: "<< dbasis.basis(0).getMaxCellLength() << "\n";
     gsInfo << "      Dofs: "<< A.numDofs() << "\n";
     // =====================================
     csvFile << "TOTAL COMPUTATIONAL TIMES, Value\n";
     csvFile << "L2 error : " << l2err << "\n";
     csvFile << "H1 error : " << h1err << "\n";
     csvFile << "H2 error : " << h2err << "\n";
     csvFile << "Mesh size: " << dbasis.basis(0).getMaxCellLength() << "\n";
     csvFile << "     Dofs: " << A.numDofs() << "\n";
     csvFile.close();
     std::cout << "Data saved to " + out + ".csv"<< std::endl; 
 
}


int main(int argc, char *argv[])
{
    real_t dt = 1e-3;
    index_t maxSteps = 10;
    index_t plotmod = 1;

    //! [Parse command line]
    bool plot = false;
    bool plot_error = false;
    index_t numRefine  = 1;
    index_t numElevate = 1;
    index_t verbose = 1;
    bool random = false;
    index_t projection_Crs = 0;
    index_t pattern = 0; // 0 for nucleation 1 for spinoidal decomposition
    std::string out("output");
    std::string fn("pde/cahn_hilliard_bvp.xml");

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addReal( "t", "dt","dt parameter",dt); // -t () or --dt ()
    cmd.addInt ( "N", "Nsteps", "Number of time steps",  maxSteps );
    cmd.addInt ( "p", "PlotMod", "Modulo for plotting",  plotmod );
    cmd.addInt ( "v", "verbose", "Verbosity level",  verbose );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("ploterror", "Create a ParaView visualization file with the projection errors", plot_error);
    cmd.addSwitch("random", "Random initial condition of the CH problem", random);
    cmd.addInt("c", "projcoars", "Projection method for coarsening", projection_Crs);
    cmd.addInt("s", "pattern", "Phase separation pattern", pattern);
    cmd.addString( "o", "output", "Output directory", out);


    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]
    gsFileData<> fd(fn);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    gsMultiPatch<> mp;
    fd.getId(0, mp); // id=0: Multipatch domain

    gsFunctionExpr<> source;
    fd.getId(1, source); // id=1: initial condition function
    // gsInfo<<"Initial condition function "<< source << "\n";

    gsBoundaryConditions<> bc;
    fd.getId(2, bc); // id=2: boundary conditions
    bc.setGeoMap(mp);
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";

    gsOptionList CHopt;
    fd.getId(3, CHopt); // id=3: reference solution

    if (random)
        gsInfo<<"Random normal initial distribution with mean "<< CHopt.askReal("mean",0.0) << " and amplitude "<< CHopt.askReal("ampl",0.005)<<"\n";
    else
        gsInfo<<"Initial condition function "<< source << "\n";

    gsOptionList TIMEopt;
    fd.getId(4, TIMEopt); // id=4: time integrator options

    gsOptionList Aopt;
    fd.getId(5, Aopt); // id=5: assembler options

    gsOptionList MESHopt;
    fd.getId(6, MESHopt); // id=6: mesher options
    //! [Read input file]

    if (mp.geoDim()==2)
        solve<2,real_t>(mp, source, bc, CHopt, TIMEopt, MESHopt, Aopt, dt, maxSteps, plotmod, plot, plot_error, numRefine, numElevate, verbose, random, projection_Crs, pattern, out);
    else if (mp.geoDim()==3)
        solve<3,real_t>(mp, source, bc, CHopt, TIMEopt, MESHopt, Aopt, dt, maxSteps, plotmod, plot, plot_error, numRefine, numElevate, verbose, random, projection_Crs, pattern, out);
    else
        GISMO_ERROR("Only 2D and 3D problems are supported.");

    return EXIT_SUCCESS;

}// end main
