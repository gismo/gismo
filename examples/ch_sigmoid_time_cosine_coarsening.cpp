/** @file cahn-hilliard.cpp

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
    ./bin/cahn-hilliard_adaptivity_example_corrected_convergence_space_hierarchical --plot -r 6 -t 0.3 -N 10 -c 0 -v 2 -s 0 -f pde/cahn_hilliard_bvp.xml 

    //  =========== CUBIC BASIS ===========
    ./bin/cahn-hilliard_adaptivity_example_corrected_convergence_space_hierarchical --plot -r 6 -t 0.3 -N 10 -c 0 -v 2 -s 0 -e 2 -f pde/cahn_hilliard_bvp.xml 
    
    //  =========== QUARTIC BASIS ===========
    ./bin/cahn-hilliard_adaptivity_example_corrected_convergence_space_hierarchical --plot -r 6 -t 0.3 -N 10 -c 0 -v 2 -s 0 -e 2 -f pde/cahn_hilliard_bvp.xml 
    
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
        gsL2Projection<real_t>::projectFunction(dbasis_fine, dbasis_coarse,fine,mp,coefs);
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

    gsStopwatch clock;

    real_t theta    = CHopt.askReal("theta",1.5);
    real_t lambda   = CHopt.askReal("lambda",1/(32*pow(EIGEN_PI,2)));
    real_t M0       = CHopt.askReal("M0",0.005);
    real_t penalty  = 1e4*lambda;

    //! [Prepare the basis]
    gsMultiBasis<> dbasis_tmp(mp,true);
    gsMultiBasis<> dbasis;

    // gsFileData<> fd1;
    // fd1.read("basis_r_7.xml");
    // fd1.getId(0,dbasis);
    // // fd1.getId(1,Coefs);

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis_tmp.setDegree( dbasis_tmp.maxCwiseDegree() + numElevate);

    if (MESHopt.askSwitch("Adaptive",true))
        MESHopt.setSwitch("THB",true);

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

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

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
    // gsFunctionExpr<> sourceQold("-9*pi*cos(6*pi*x)*cos(6*pi*y)*((2*sin((2*pi*t)/3))/27 - (4930115259914363*pi*cos((2*pi*t)/3))/8796093022208 + pi*cos((2*pi*t)/3)^3*(24*cos(6*pi*y)^2 - cos(6*pi*x)^2*(72*cos(6*pi*y)^2 - 24)))",2);
    // gsFunctionExpr<> sourceQnew("-9*pi*cos(6*pi*x)*cos(6*pi*y)*((2*sin((2*pi*t)/3))/27 - (4930115259914363*pi*cos((2*pi*t)/3))/8796093022208 + pi*cos((2*pi*t)/3)^3*(24*cos(6*pi*y)^2 - cos(6*pi*x)^2*(72*cos(6*pi*y)^2 - 24)))",2);
    // Source (volume integral) function for manufactured solution cos(pi*x) * cos(pi*y) * cos(2/3*pi*t)
    gsFunctionExpr<> sourceQold("-9*pi*cos(pi*x)*cos(pi*y)*((2*sin((2*pi*t)/3))/27 - (274134357077347*pi*cos((2*pi*t)/3))/1266637395197952 + pi*cos((2*pi*t)/3)^3*((2*cos(pi*y)^2)/3 - cos(pi*x)^2*(2*cos(pi*y)^2 - 2/3)))",2);
    gsFunctionExpr<> sourceQnew("-9*pi*cos(pi*x)*cos(pi*y)*((2*sin((2*pi*t)/3))/27 - (274134357077347*pi*cos((2*pi*t)/3))/1266637395197952 + pi*cos((2*pi*t)/3)^3*((2*cos(pi*y)^2)/3 - cos(pi*x)^2*(2*cos(pi*y)^2 - 2/3)))",2);


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
    // gsFunctionExpr<> bc1("-18*pi*cos((2*pi*t)/3)*cos(6*pi*y)*sin(6*pi*x)*(cos((2*pi*t)/3)^2*cos(6*pi*x)^2*cos(6*pi*y)^2 + 927703176743281/42221246506598400)","-18*pi*cos((2*pi*t)/3)*cos(6*pi*x)*sin(6*pi*y)*(cos((2*pi*t)/3)^2*cos(6*pi*x)^2*cos(6*pi*y)^2 + 927703176743281/42221246506598400)",2);
    
    // For manufactured solution cos(pi*x) * cos(pi*y) * cos(2/3*pi*t)
    gsFunctionExpr<> bc1("-3*pi*cos((2*pi*t)/3)*cos(pi*y)*sin(pi*x)*(cos((2*pi*t)/3)^2*cos(pi*x)^2*cos(pi*y)^2 + 274134357077347/844424930131968)", "-3*pi*cos((2*pi*t)/3)*cos(pi*x)*sin(pi*y)*(cos((2*pi*t)/3)^2*cos(pi*x)^2*cos(pi*y)^2 + 274134357077347/844424930131968)",2);

    // (2) Laplace boundary condition
    // gsFunctionExpr<> bc2("-72*pi^2*cos((2*pi*t)/3)*cos(6*pi*x)*cos(6*pi*y)",2); // should be correct
    
    // For manufactured solution cos(pi*x) * cos(pi*y) * cos(2/3*pi*t)
    gsFunctionExpr<> bc2("-3*pi*cos((2*pi*t)/3)*cos(pi*y)*sin(pi*x)*(cos((2*pi*t)/3)^2*cos(pi*x)^2*cos(pi*y)^2 + 274134357077347/844424930131968)", "-3*pi*cos((2*pi*t)/3)*cos(pi*x)*sin(pi*y)*(cos((2*pi*t)/3)^2*cos(pi*x)^2*cos(pi*y)^2 + 274134357077347/844424930131968)",2);

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
    }
    else
    {
        // %%%%%%%%%%%%%%%%%%%%%%%% Analytical intial condition %%%%%%%%%%%%%%%%%%%%%%%%
        GISMO_ASSERT(mp.geoDim()==source.domainDim(),"Domain dimension of the source function should be equal to the geometry dimension, but "<<source.domainDim()<<"!="<<mp.geoDim());
        gsMatrix<> tmp;
        Cold.setZero(A.numDofs(),1);
        gsFunctionExpr<> initial_cond("tanh(20*x + 2*cos((pi*t)/3) - 10)",2);
        real_t error = gsL2Projection<real_t>::projectFunction(dbasis,initial_cond,mp,tmp);  // 3rd arg has to be multipatch
        if (verbose>0) gsInfo << "L2 projection error "<<error<<"\n";
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
    csvFile.open(out+"/dofs.csv");
    csvFile << "TimeStep, NumDOFs, Mass, ErrorRefCnew, ErrorRefdCnew, ErrorCrsCnew, ErrorCrsdCnew\n";

    gsVector<> pt(2,1); pt<<0.5, 0.5;

    // auto u_manufactured = ev.getVariable(source_time, G);
    // gsInfo<< "Value at mid-point :" << ev.eval(neumann_manuf, pt) <<"\n";

    // real_t mass2 = ev.integral(meas(G)*cold);
    // gsInfo<<mass2<<"\n";
    real_t told, tnew;

    
    // =============================================================================================
    // Do the coarsening based on the analytical solution!! (initial condition)
    
    if (MESHopt.askSwitch("Adaptive",true))
        mesher.rebuild();

    
    for (index_t step = 0; step!=maxSteps; step++)
    {
        gsInfo<< "Coarsening step: "<<step<<"\n";
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

        gsDebugVar(coarsen.totalSize());

        // If elements are marked for refinement
        if (coarsen.totalSize()!=0)
        {
            // Refine dbasis
            if (verbose>1) gsInfo<<"Basis before coarsening:\n "<<dbasis.basis(0)<<"\n";
            mesher.unrefine(coarsen);
            if (verbose>1) gsInfo<<"Basis after coarsening:\n "<<dbasis.basis(0)<<"\n";
        }// coarsen
    }


    if (plot)
    {
        // Export the mesh
        collection.newTimeStep(&mp);
        // collection.addField(cnew,"numerical solution");
        gsFunctionExpr<> source_time("tanh(20*x + 2*cos((pi*t)/3) - 10)",2);
        source_time.set_t(0);
        auto u_manufactured = ev.getVariable(source_time, G);
        collection.addField(u_manufactured,"analytical solution");
        // collection.addField((u_manufactured - cnew_sol).sqNorm(), "L2 error");
        // collection.addField(u_manufactured.sqNorm() - cnew_sol.sqNorm(), "difference");
        collection.saveTimeStep();
    }

    if (plot)
        collection.save(); 

    // =================== Save the coarsened basis !!! ===================
    gsFileData<> fdbas;
    fdbas.add(dbasis,0);
    fdbas.save("basis_sigmoid_mark_"+std::to_string(MESHopt.getReal("CoarsenParam"))+"_r_"+std::to_string(numRefine)+".xml");
    gsInfo << "Exported to "<<  "basis_sigmoid_mark_"+std::to_string(MESHopt.getReal("CoarsenParam"))+"_r_"+std::to_string(numRefine)+".xml" <<"\n";
    A.initSystem();
    gsInfo << "System size: "<< A.numDofs() <<"\n";
    gsInfo << "Basis size: "<< dbasis.size() <<"\n";
    // ====================================================================

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
