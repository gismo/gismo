/** @file cahn-hilliard.cpp

    @brief Tutorial on how to use expression assembler to solve the Cahn-Hilliard equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Marsala (UniFi)
               H.M. Verhelst (UniFi)
               L. Venta Viñuela (UniPv)


    2D
    ./build_pardiso/bin/cahn-hilliard_adaptivity_paper_examplemark_plotmesh -r 6 -t 1.5e-3 -N 1000 -c 0 -v 2 -s 0 -p 20 --plot           --random -f pde/cahn_hilliard_bvp_Sameelement1.xml      -o “output/TP_nucleation_2D”
    ./build_pardiso/bin/cahn-hilliard_adaptivity_paper_examplemark_plotmesh -r 6 -t 1.5e-3 -N 1000 -c 1 -v 2 -s 0 -p 20 --plot --plotmesh --random -f pde/cahn_hilliard_bvp_ad_Sameelement1.xml  -o “output/QI_nucleation_2D”
    3D
    ./build_pardiso/bin/cahn-hilliard_adaptivity_paper_examplemark_plotmesh  -r 6 -t 1.5e-3 -N 1000 -c 0 -v 2 -s 0 -p 20 --plot           --random -f pde/cahn_hilliard_3d_bvp.xml      -o "output/TP_nucleation_3D"
    ./build_pardiso/bin/cahn-hilliard_adaptivity_paper_examplemark_plotmesh  -r 6 -t 1.5e-3 -N 1000 -c 1 -v 2 -s 0 -p 20 --plot --plotmesh --random -f pde/cahn_hilliard_3d_bvp_ad.xml  -o “output/QI_nucleation_3D”

    -----------------------------------------------------------------------
    TODO;
    - Change hmax to a gsExprAssembler<>::element el; el.diam();
    -----------------------------------------------------------------------

*/

//! [Include namespace]
#include <gismo.h>
#include <iomanip> // for std::setw, std::setfill, std::setprecision
#include <gsHSplines/gsHElementMarker.h>

using namespace gismo;
//! [Include namespace]

template <class T>
T computeNorm(const gsMultiBasis<T>  & integrationBasis,
                const gsFunctionSet<T> & currentBasis,
                const gsFunctionSet<T> & geometryMap,
                const gsFunctionSet<T> & sourceFunction1,
                const gsMatrix<T> & coefs,
                const gsOptionList     & options)
{
  
    // Create an assembler
    gsExprAssembler<T> A(1,1);
    A.options().update(options,gsOptionList::addIfUnknown); 
    // gsInfo<< "Assembler options: " << A.options() << "\n";

    // Set the integration elements
    A.setIntegrationElements(integrationBasis);

    auto u = A.getSpace(currentBasis);
    u.setup();
    // gsMatrix<T> tmp = coefs; // coefs is const
    // auto sol= A.getSolution(u,tmp);

    A.initSystem();

    // Assign the geometry map
    typename gsExprAssembler<T>::geometryMap G = A.getMap(geometryMap);

    gsExprEvaluator<T> ev(A);
    // ev.options().update(A.options(),gsOptionList::addIfUnknown);

    auto sol_before = ev.getVariable(sourceFunction1,G); // solution before projection      
    gsGeometry<>::uPtr sol_after_ptr  = currentBasis.basis(0).makeGeometry(coefs);
    auto sol_after = ev.getVariable(*sol_after_ptr,G);

    // gsInfo<<currentBasis.basis(0)<<"\n";
    // ev.options().setSwitch("SameElement",0);        

    return ev.integral((sol_before-sol_after).sqNorm() * meas(G));

}

template <class T, class Expr>
T computeInt(const gsMultiBasis<T>  & integrationBasis,
                const gsFunctionSet<T> & geometryMap,
                const Expr & sourceFunction1,
                const gsOptionList     & options)
{
  
    // Create an assembler
    gsExprAssembler<T> A(1,1);
    A.options().update(options,gsOptionList::addIfUnknown); 

    // Set the integration elements
    A.setIntegrationElements(integrationBasis);

    // Assign the geometry map
    typename gsExprAssembler<T>::geometryMap G = A.getMap(geometryMap);

    gsExprEvaluator<T> ev(A);
    ev.options().update(A.options(),gsOptionList::addIfUnknown);
    ev.options().setSwitch("SameElement",false); // add warning to remind that im deactivating sameelemnt        

    return ev.integral(sourceFunction1 * meas(G));

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
            bool & plot_mesh,
            index_t & numRefine,
            index_t & numElevate,
            index_t & verbose,
            bool & random,
            index_t & projection_Crs,
            index_t & pattern,
            index_t & timemethod,
            std::string out)
{
    typedef typename gsHElementHelper<dim,real_t>::element_t element_t;
    typedef typename gsHElementHelper<dim,real_t>::HElementContainer HElementContainer;

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
    index_t maxRefIt;


    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    A.options().setSwitch("SameElement",Aopt.askSwitch("SameElement",true));
#ifdef GISMO_WITH_PARDISO
    A.options().addString("LinearSolver", "Name of the linear solver to be used", "PardisoLU");
#   else
    A.options().addString("LinearSolver", "Name of the linear solver to be used", "SimplicialLDLT");
#endif
    A.options().setReal("quA", Aopt.askReal("quA",1));
    A.options().setInt("quB", Aopt.askInt("quB",1));


    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    // A.setIntegrationElements(dbasis);
    // =====================================================================

    gsMultiBasis<> ibasis, evbasis;
    // Cast every basis of dbasis to a gsTHBSplineBasis
    for (size_t p=0; p!=dbasis.nBases(); p++)
    {
        GISMO_ASSERT((dynamic_cast<gsTHBSplineBasis<dim,real_t>*>(&dbasis.basis(p))), "Basis is not a gsTHBSplineBasis");
        gsTHBSplineBasis<dim,real_t> & b = static_cast<gsTHBSplineBasis<dim,real_t> & >(dbasis.basis(p));
        ibasis.addBasis(b.tensorLevel(b.maxLevel()).clone());
    }
    for (size_t p=0; p!=dbasis.nBases(); p++)
    {
        GISMO_ASSERT((dynamic_cast<gsTHBSplineBasis<dim,real_t>*>(&dbasis.basis(p))), "Basis is not a gsTHBSplineBasis");
        gsTHBSplineBasis<dim,real_t> & b = static_cast<gsTHBSplineBasis<dim,real_t> & >(dbasis.basis(p));
        evbasis.addBasis(b.tensorLevel(b.maxLevel()).clone());
    }
    A.setIntegrationElements(ibasis); // USE FINEST TENSOR BASIS FOR INTEGRATION
    gsExprEvaluator<> ev(A);
    ev.options().update(A.options(),gsOptionList::addIfUnknown); //?? do I have to do this?
    A.options().update(ev.options(),gsOptionList::addIfUnknown); 
    ev.options().setSwitch("SameElement",false);
    gsDebugVar(dbasis.basis(0).maxDegree());
    gsDebugVar(ibasis.basis(0));

    // Set the geometry map
    // geometryMap G = A.getMap(surface);
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
    // auto     cnew  = A.getCoeff(mp_cold );
    // auto     dcnew = A.getCoeff(mp_dcold);

    // solution cold = A.getSolution(w, Cold); // Cold
    solution cnew_sol = A.getSolution(w, Cnew); // Cnew
    // solution dcold = A.getSolution(w, dCold); // \dot{Cold}
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

    // // Solution variables for the intermediate solutions (during time integration)
    // solution c = A.getSolution(w, Calpha); // C
    // solution dc = A.getSolution(w, dCalpha); // \dot{C}


    // Derivatives of the polynomial double well potential (M. Kästner et al., 2016)
    auto dmu_c = - 1.0 + 3.0 * (c*c).val(); // f_2 (second derivative of double well)
    auto ddmu_c = 6*c.val(); // f_3 (third derivative of double well)

    // Mobility
    auto M_c  = 1.0 + 0.0*c.val(); // replace with const_expr(1.0) instead of using 0*c
    auto dM_c = 0.0 * gradc; // replace with const_expr(1.0) instead of using 0*c!!

    auto residual = w*dc + // M
                    M_c.val() * igrad(w,G) * dmu_c * gradc.tr() + // F_bar
                    M_c.val() * ilapl(w,G) *lambda *laplc.val(); // K_laplacian
                    // lambda*laplc.val()*igrad(w,G)*dM_c.tr() + // term gradient mobility!

    //! [Problem setup]

    // ![Initialize the assembler]
    w.setup(bc, dirichlet::l2Projection, 0);
    A.initSystem();
    // ![Initialize the assembler]
gsSparseSolver<>::uPtr solver;
#ifdef GISMO_WITH_PARDISO
    solver = gsSparseSolver<>::get("PardisoLU");
#   else
    solver = gsSparseSolver<>::get("LU");
#endif

    gsMatrix<> Q;
    gsSparseMatrix<> K, K_const, K_mass, K_mass_full;

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
        // gsMatrix<> tmp = gsMatrix<>::Random(dbasis.size(),1);
        // tmp.array() *= CHopt.askReal("ampl",0.005); //random uniform variable in [-0.05,0.05]
        // tmp.array() += CHopt.askReal("mean",0.0); // 0.45
        // gsDebugVar(tmp);
        // mp_cold.addPatch(dbasis.basis(0).makeGeometry(tmp));
        // tmp.setZero();
        // mp_dcold.addPatch(dbasis.basis(0).makeGeometry(tmp));
        // gsInfo<<mp_cold.patch(0).coefs().size()<<"\n";
        // gsInfo<<mp_dcold.patch(0).coefs().size()<<"\n";



        // gsFileData<> fd2;
        // fd2.add(mp_cold,0); // mp
        // fd2.save("multipatch_3d.xml");
        // gsInfo << "Exported to multipatch_3d.xml \n";
        

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
    
        gsFileData<> fd1;
        std::string file_name;
        if (pattern==0) // nucleation
        {
            if (dim==2)
                file_name = "pde/ic_nucleation.xml";
            else if (dim==3)
                file_name = "pde/ic_nucleation_3d.xml";
        }
        else if (pattern==1) // spinoidal
        {
            if (dim==2)
                file_name = "pde/ic_spinoidal.xml";
            else if (dim==3)
                file_name = "pde/ic_spinoidal_3d.xml";
        }

        gsMultiPatch<> MP;
        fd1.read(file_name);
        fd1.getId(0,MP);

        mp_cold.addPatch(MP.patch(0));
        Cold.setZero(dbasis.size(),1);
        mp_dcold.addPatch(dbasis.basis(0).makeGeometry(Cold));
        gsInfo<<"Imported multipatch\n";
    }
    else
    {
        // %%%%%%%%%%%%%%%%%%%%%%%% Analytical intial condition %%%%%%%%%%%%%%%%%%%%%%%%
        GISMO_ASSERT(mp.geoDim()==source.domainDim(),"Domain dimension of the source function should be equal to the geometry dimension, but "<<source.domainDim()<<"!="<<mp.geoDim());
        gsMatrix<> tmp;
        Cold.setZero(A.numDofs(),1);
        real_t error = gsL2Projection<real_t>::project(dbasis,source,mp,tmp);  // 3rd arg has to be multipatch
        if (verbose>0) gsInfo << "L2 projection error "<<error<<"\n";
        mp_cold.addPatch(dbasis.basis(0).makeGeometry(tmp));
        tmp.setZero();
        mp_dcold.addPatch(dbasis.basis(0).makeGeometry(tmp));
    }

    mp_cnew = mp_cold;
    mp_dcnew = mp_dcold;

    gsInfo<< mp_cold.patch(0).coefs().maxCoeff() << "\n";
    gsInfo<< mp_cold.patch(0).coefs().minCoeff() << "\n";
    gsInfo<< mp_dcold.patch(0).coefs().maxCoeff()  << "\n";
    gsInfo<< mp_dcold.patch(0).coefs().minCoeff()  << "\n";

    real_t Q0norm = 1, Qnorm = 10;
    real_t tol = TIMEopt.askReal("tol",1e-4);

    gsFileManager::mkdir(out); // create output folder
    std::string path_mesh = out + "/mesh_pvd"; // create output folder for mesh
    gsInfo << "Mesh directory path: " << path_mesh << "\n";
    GISMO_ENSURE(gsFileManager::mkdir(path_mesh),"Error making directory "+path_mesh);

    gsParaviewCollection meshcollection(out+"/mesh");
    int plotID = 0;

    gsParaviewCollection collection(out+"/solution", &ev);
    collection.options().setSwitch("plotElements", false);
    collection.options().setInt("plotElements.resolution", 4); //amount of elements (segment?)
    collection.options().setInt("numPoints",(mp.geoDim()==3) ? 100000 : 5000);

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
    real_t error_proj_c = 0;
    real_t error_proj_dc = 0;

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

    gsHElementMarker<dim,real_t> marker(dbasis.basis(0));
    marker.options().setInt("MaxLevel",numRefine);
    marker.options().setInt("RefineRule",1);
    marker.options().setReal("RefineParam",0.1); // is this correct?
    marker.options().setInt("CoarsenRule",1);
    marker.options().setReal("CoarsenParam",0.1); // is this correct?
    marker.options().setSwitch("Absolute",true);

    // ! [Load mesher options]

    std::ofstream csvFile;
    csvFile.open(out+"/output_data.csv");
    csvFile << "TimeStep,"
            << "RefIt,"
            << "NumDOFs,"
            << "nSolvesStep,"
            << "Mass,"
            << "Bulk_energy,"
            << "Int_energy,"
            << "Total_energy,"
            << "proj_err c,"
            << "proj_err dc,"
            << "AT,"
            << "ST,"
            << "PT,"
            << "MarkRefT,"
            << "RefT,"
            << "MarkCrsT,"
            << "CrsT\n";

    std::ofstream csvFile2;
    csvFile2.open(out+"/elements.csv");
    csvFile2 << std::left
            << "TimeStep,"          // sin std::setw: un CSV no necesita espacios
            << "RefIt,"             
            << "NumDOFs,"           
            << "nSolvesStep";

    std::vector<size_t> nElements(numRefine+1,0); // size equal to the number of lebels    

    for (index_t i = 0; i < nElements.size(); ++i) {
        csvFile2 << ",Elements (lvl " << i << ")";
    }
    
    csvFile2 <<","<< "MarkedRef"; // Marked elements for refinement
    csvFile2 <<","<< "tot_elems"; // 
    csvFile2 <<","<< "MarkedRef/tot_elems"; // 
    csvFile2 << '\n';   

    gsVector<> pt(2,1); pt<<0.5, 0.5;

    real_t assemblyTimestep, solverTimestep, projTimestep;
    real_t assemblyTimeRefIt, solverTimeRefIt, projTimeRefIt;
    real_t markRefTimeRefIt, refTimeRefIt, markCrsTimeStep, crsTimeStep;  
    real_t markRefTimeStep, refTimeStep;

    real_t crsTime = 0;
    real_t markCrsTime = 0;
    real_t refTime = 0;
    real_t markRefTime = 0;

    std::string method;

    // Set parameters for time integration method
    if (timemethod == 0)
    {
        tmp_alpha_m = tmp_alpha_f = tmp_gamma = 1;
        tmp_alpha_m_func.setValue(tmp_alpha_m,dim);
        tmp_alpha_f_func.setValue(tmp_alpha_f,dim);
        method = "Backward Euler ";
    }
    else if (timemethod == 1) // Generalized-alpha method
    {
        tmp_alpha_m = alpha_m;
        tmp_alpha_f = alpha_f;
        tmp_gamma   = gamma;
        tmp_alpha_m_func.setValue(tmp_alpha_m,dim);
        tmp_alpha_f_func.setValue(tmp_alpha_f,dim);
        method = "Generalized Alpha ";
    }

    // We do this once to be able to get the coefficients of the initial condition from the xml file
    if (MESHopt.askSwitch("Adaptive",true)==0)
    {
        // Tensor product solution (to be able to have finer meshes than the one from the xml file)
        gsL2Projection<real_t>::project(dbasis,ibasis,mp,mp_cold.patch(0),CnewF,A.options());
        gsL2Projection<real_t>::project(dbasis,ibasis,mp,mp_dcold.patch(0),dCnewF,A.options());
    }

    for (index_t step = 0; step!=maxSteps; step++)
    {
        assemblyTimestep = 0;
        solverTimestep   = 0;
        projTimestep     = 0;

        markRefTimeStep = 0;
        refTimeStep     = 0;
        markCrsTimeStep = 0;
        crsTimeStep     = 0;
 
        // if ((step >= 80 && step <= 200) || (step >= 320 && step <= 460) || (step >= 700 && step <= 800))
            maxRefIt = MESHopt.askInt("RefIt",5);
        // else
        //     maxRefIt = 1;

        for (index_t refIt = 0; refIt!=maxRefIt; refIt++)
        {
            
            nSolvesStep       = 0;
            assemblyTimeRefIt = 0;
            solverTimeRefIt   = 0;
            projTimeRefIt     = 0;
            
            markRefTimeRefIt  = 0;
            refTimeRefIt      = 0;

            // ===== Compute % elements in each level of the hierarchy =====
            GISMO_ASSERT(dbasis.nBases()==1, "Implemented for only one basis");
            // Assigns errors to the elements, based on their level: e = 10^{l-1}
            typename gsBasis<>::domainIter domIt = dbasis.basis(0).domain()->beginAll();
            typename gsBasis<>::domainIter domItEnd = dbasis.basis(0).domain()->endAll();
            gsHDomainIterator<real_t,dim> * domHIt = dynamic_cast<gsHDomainIterator<real_t,dim> *>(domIt.get());
            std::vector<size_t> nElements(numRefine+1,0); // size equal to the number of lebels    
            std::vector<real_t> elAreas(numRefine+1,0); // size equal to the number of lebels
            for (; domIt<domItEnd; ++domIt)
            {
                nElements[domHIt->getLevel()]++; 
                // elAreas[domHIt->getLevel()] += (domHIt->upperCorner()-domHIt->lowerCorner()).prod();
                elAreas[domHIt->getLevel()] += domHIt->volume();
            }
            GISMO_ENSURE(math::abs(gsAsVector<real_t>(elAreas).sum()-1.0) < 1e-12, "Element area should be close to 1.0, but is " << gsAsVector<real_t>(elAreas).sum());

            csvFile2 << step         << ','  
                    << refIt         << ','   
                    << A.numDofs()   << ','   
                    << nSolvesStep;         
            for (index_t i = 0; i < elAreas.size(); ++i)
                csvFile2 << ',' << elAreas[i];
            if (MESHopt.askSwitch("Adaptive",true)==0)
            {
              csvFile2 << '\n';                
              csvFile2.flush();     
            }
            // ================================================================================

            clock.restart();
            if (MESHopt.askSwitch("Adaptive",true))
            {
                if (projection_Crs == 0)
                {
                    // Add if for mp_cnew!!!
                    // error_proj_c = gsL2Projection<real_t>::project(dbasis,ibasis,mp,mp_cold.patch(0),CnewF,A.options());
                    // error_proj_dc = gsL2Projection<real_t>::project(dbasis,ibasis,mp,mp_dcold.patch(0),dCnewF,A.options());    
                    // gsDebug<<"Error in the L2 projection of the initial condition (c) : "<<error_proj_c<<"\n";
                    // gsDebug<<"Error in the L2 projection of the initial condition (dc): "<<error_proj_dc<<"\n";  

                    gsL2Projection<real_t>::project(dbasis,ibasis,mp,mp_cold.patch(0),CnewF,A.options());
                    gsL2Projection<real_t>::project(dbasis,ibasis,mp,mp_dcold.patch(0),dCnewF,A.options()); 

                    // l2_proj(dbasis,ibasis,mp,mp_cold.patch(0),K_mass,CnewF,A.options());
                    // l2_proj(dbasis,ibasis,mp,mp_dcold.patch(0),K_mass,dCnewF,A.options()); 
                                    
                    // // Projection of solution c
                    // auto f1 = A_proj.getCoeff(mp_cold.patch(0));
                    // A_proj.initSystem();
                    // A_proj.assemble(u * f1 * meas(G)); // assemble rhs
                    // solver->compute(K_mass);
                    // CnewF = solver->solve(A_proj.rhs());

                    // // Projection of velocity field dc
                    // auto f2 = A_proj.getCoeff(mp_dcold.patch(0));
                    // A_proj.initSystem();
                    // A_proj.assemble(u * f2 * meas(G)); // assemble rhs
                    // solver->compute(K_mass);
                    // dCnewF = solver->solve(A_proj.rhs());
                }
                else if (projection_Crs == 1)
                {
                    gsQuasiInterpolate<real_t>::localIntpl(dbasis.basis(0),mp_cold.patch(0),CnewF);
                    gsQuasiInterpolate<real_t>::localIntpl(dbasis.basis(0),mp_dcold.patch(0),dCnewF);    
                }
                else if (projection_Crs == 2)
                {
                    gsQuasiInterpolate<real_t>::Schoenberg(dbasis.basis(0),mp_cold.patch(0),CnewF);
                    gsQuasiInterpolate<real_t>::Schoenberg(dbasis.basis(0),mp_dcold.patch(0),dCnewF);    
                }
            }
            else if (step > 0)
            {
                // Tensor product solution (to be able to have finer meshes than the one from the xml file)
                // gsL2Projection<real_t>::project(dbasis,ibasis,mp,mp_cold.patch(0),CnewF,A.options());
                // gsL2Projection<real_t>::project(dbasis,ibasis,mp,mp_dcold.patch(0),dCnewF,A.options());
                // Tensor product (mesh stays the same -- no projection needed)
                CnewF = mp_cold.patch(0).coefs();
                dCnewF = mp_dcold.patch(0).coefs();
            }
            projTimeRefIt += clock.stop();

            real_t error_proj_c  = computeNorm(evbasis, dbasis, mp, mp_cold.patch(0), CnewF, A.options());
            real_t error_proj_dc = computeNorm(evbasis, dbasis, mp, mp_dcold.patch(0), dCnewF, A.options());

            if (step>=128)
            {
                gsInfo<< "error_proj_c: " << error_proj_c << "\n";
                gsInfo<< "error_proj_dc: " << error_proj_dc << "\n";
                gsInfo<< "Norm cold: " << CnewF.norm() << "\n";
                gsInfo<< "Norm cold: " << dCnewF.norm() << "\n";
                gsInfo << "Cnew head: " << Cnew.topRows(10).transpose() << "\n";
                gsInfo << "dCnew head: " << dCnew.topRows(10).transpose() << "\n";
                gsInfo << "w.mapper() mapSize: " << w.mapper().mapSize() << "\n";
                gsInfo << "w.mapper() freeSize: " << w.mapper().freeSize() << "\n";
            }
                      // // =========== Compute error of the projection of both variables ===========
            // if (step >= 102)
            // { 
            //     real_t error_proj_c  = computeNorm(ibasis, dbasis, mp, mp_cold.patch(0), CnewF, A.options());
            //     real_t error_proj_dc = computeNorm(ibasis, dbasis, mp, mp_dcold.patch(0), dCnewF, A.options());
            //     gsInfo << "=====================================================================\n";
            //     gsInfo << "Step: " << step << ", RefIt: " << refIt << "\n";
            //     gsInfo << "Projection error of c: " << error_proj_c << "\n";
            //     gsInfo << "Projection error of dc: " << error_proj_dc << "\n";
            //     gsInfo << "Basis mp_cold size: " << mp_cold.basis(0).size() << "\n";
            //     gsInfo << "Basis dbasis size: " << dbasis.basis(0).size() << "\n";
            //     gsInfo << "Basis mp_cold: "<< mp_cold.basis(0) << "\n";
            //     gsInfo << "Basis dbasis: "<< dbasis.basis(0) << "\n";
            //     gsInfo << "Max CnewF: " << CnewF.maxCoeff() << "\n";
            //     gsInfo << "Min CnewF: " << CnewF.minCoeff() << "\n";
            //     gsInfo << "Max dCnewF: " << dCnewF.maxCoeff() << "\n";
            //     gsInfo << "Min dCnewF: " << dCnewF.minCoeff() << "\n";
            //     gsInfo << "=====================================================================\n";
            // }
            // // ==========================================================================
            
            // =========================Assemble constant terms===============================
            // Resize the data structure inside the mesher
            if (MESHopt.askSwitch("Adaptive",true))
                mesher.rebuild();

            // Setup the space (compute Dirichlet BCs)
            w.setup(bc, dirichlet::l2Projection, 0);

            // Reset the assembler
            A.initSystem();
            //  Assemble all the terms that are independent of the solution !!!! Mass matrix and bilaplacian
            clock.restart();
            A.assemble(meas(G) * (w*w.tr()*tmp_alpha_m +// K_m
                                (tmp_alpha_f * tmp_gamma * dt)* (lambda * ilapl(w,G) * ilapl(w,G).tr()))); // K_laplacian
            // A.assemble(meas(G) * ((tmp_alpha_f * tmp_gamma * dt)* (lambda * ilapl(w,G) * ilapl(w,G).tr()))); // K_laplacian

            K_const = A.giveMatrix();
            // K_const = tmp_alpha_m*K_mass + A.giveMatrix(); // .giveMatrix() moves the matrix A into K_const (avoids having two matrices A and K_const)
            assemblyTimeRefIt += clock.stop();
            // ===============================================================================

            gsInfo<< "Number of Elements: " << dbasis.basis(0).numElements() << "\n";
            
        
            // if (step==203)
            // {
            //     gsInfo<< "ibasis size: " << ibasis.basis(0).size() << "\n";
            //     gsInfo<< "dbasis size: " << dbasis.basis(0).size() << "\n";
            //     gsInfo<< "ibasis size: " << ibasis.basis(0).numElements() << "\n";
            //     gsInfo<< "dbasis size: " << dbasis.basis(0).numElements()  << "\n";
            //     gsMesh<> mesh0(dbasis.basis(0));
            //     gsMesh<> mesh1(ibasis.basis(0));
            //     gsInfo<<"Exporting mesh to" <<out+"/trial_dbasis_mesh_"+std::to_string(step)<<"\n";
            //     gsWriteParaview(mesh0, out+"/trial_dbasis_mesh_"+std::to_string(step));
            //     gsWriteParaview(mesh1, out+"/trial_ibasis_mesh_"+std::to_string(step));
            // }

            for (index_t dt_it = 0; dt_it != lmax; dt_it++)
            {
                if (verbose>0) gsInfo<<"Time step "<<step<<"/"<<maxSteps<<", iteration "<<dt_it<<": dt = "<<dt<<", [t_start,t_end] = ["<<time<<","<<time+dt<<"]"<<"\n";
                converged = false;


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


                Q0norm = 1;
                Qnorm = 10;

                // ==== For Nitsche BC ====
                Cold.setZero(A.numDofs(),1);
                dCold.setZero(A.numDofs(),1);
                Cold = Cnew;
                dCold = dCnew;
                // ========================

                for (index_t it = 0; it!= maxIt; it++)
                {
                    // if (step>=102)
                    // {
                    //     gsExprEvaluator<> ev1(A);
                    //     ev1.options().update(A.options(),gsOptionList::addIfUnknown); //?? do I have to do this?
                    //     ev1.options().setSwitch("SameElement",false);

                    //     gsInfo << "=====================================================================\n";
                    //     A.initSystem();
                    //     A.assemble(w*(cold));
                    //     gsInfo<< ev1.options() <<"\n";

                    //     gsInfo<< "Norm of cold: "<<A.rhs().norm()<<"\n";
                    //     // A.initSystem();
                    //     // A.assemble(w*(cnew_sol));
                    //     // gsInfo<< "Norm of cnew_sol: "<<A.rhs().norm()<<"\n";
                    //     A.initSystem();
                    //     A.assemble(w*(dcold));
                    //     gsInfo<< "Norm of dcold: "<<A.rhs().norm()<<"\n";
                    //     gsInfo<< "Int(cold): "<< ev1.integral(cold) <<"\n";
                    //     gsInfo<< "Int(dcold): "<< ev1.integral(dcold) <<"\n";
                    //     gsInfo<< "Domain elements A "<<  A.domain().numElements() <<"\n";
                    //     gsInfo<< "dbasis elements " << dbasis.basis(0).numElements()<<"\n";
                    //     gsInfo<< "ibasis elements " << ibasis.basis(0).numElements()<<"\n";
                    //     // A.initSystem();
                    //     // A.assemble(w*(dcnew_sol));
                    //     // gsInfo<< "Norm of dcnew_sol: "<<A.rhs().norm()<<"\n";
                        
                    //     // A.initSystem();
                    //     // A.assemble(w*(dcnew_sol-dcold));
                    //     // gsInfo<< "Norm of dcnew_sol-dcold: "<<A.rhs().norm()<<"\n";
                        
                    //     // A.initSystem();
                    //     // A.assemble(w*(cnew_sol-cold));
                    //     // gsInfo<< "Norm of cnew_sol-cold: "<<A.rhs().norm()<<"\n";
                    //     gsInfo << "=====================================================================\n";

                    // }
                    
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
                    assemblyTimeRefIt += clock.stop();

                    Q = A.rhs();

                    // Assemble the Nitsche BC on the sides with Neumann condition
                    // A.clearMatrix(); // Resets to zero the values of the already allocated to matrix (LHS)
                    A.initMatrix();
                    clock.restart();
                    A.assembleBdr(bc.get("Neumann"), - lambda * igrad(w,G) *  nv(G)  * ilapl(w,G).tr() + // consistency term
                                penalty * (igrad(w,G) * nv(G).normalized()) * hmax * (igrad(w,G) * nv(G)).tr() - // penalty (stabilizing) term
                                lambda * ilapl(w,G) * (igrad(w,G)  * nv(G)).tr()); // symmetry term
                    assemblyTimeRefIt += clock.stop();

                    K_nitsche = A.giveMatrix(); // .giveMatrix() moves the matrix A into K_nitche (avoids having two matrices A and K_nitsche)

                    if (bc.get("Neumann").size()!=0)
                        Q.noalias() += K_nitsche * Calpha; // add the residual term from Nitche (using the matrix )

                    // Check the convergence conditions
                    if (it == 0) Q0norm = Q.norm();
                    else         Qnorm = Q.norm();

                    if (verbose==2) gsInfo<<"\t\tNR iter   "<<it<<": res = "<<Qnorm/Q0norm<<"\n";

                    if (it>0 && Qnorm/Q0norm < tol)
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

                    A.initMatrix();
                    // Assembly of the tangent stiffness matrix (K_m and K_f simultaneously) %%
                    clock.restart();
                    // A.assemble(meas(G) * (w*w.tr()*tmp_alpha_m +// K_m
                    //                     (tmp_alpha_f * tmp_gamma * dt)* (dmu_c *igrad(w,G) * igrad(w,G).tr() + // K_f1
                    //                     ddmu_c * igrad(w,G) * gradc.tr() * w.tr() + // K_f2
                    //                     lambda * ilapl(w,G) * ilapl(w,G).tr()))); // K_laplacian
                    //                     // lambda * igrad(w,G)*dM_c.tr()*ilapl(w,G).tr()   +  // K_mobility
                    A.assemble(meas(G) * (tmp_alpha_f * tmp_gamma * dt)* (dmu_c *igrad(w,G) * igrad(w,G).tr() + // K_f1
                                        ddmu_c * igrad(w,G) * gradc.tr() * w.tr())); // K_laplacian
                    assemblyTimeRefIt += clock.stop();

                    K = A.giveMatrix();
                    K += K_const; // add the constant part of the stiffness matrix (K_m and K_laplacian)

                    if (bc.get("Neumann").size()!=0)
                        K += (tmp_alpha_f * tmp_gamma * dt) * K_nitsche; // add the Nitsche term to the stiffness matrix


                    clock.restart();
                    
                    solver->compute(K); // K + K_linear ?? 

                    dCupdate = solver->solve(-Q);
                    solverTimeRefIt += clock.stop();
                    nSolves++;
                    nSolvesStep++;

                    dCnew += dCupdate;
                    Cnew.noalias() += (tmp_gamma*dt)*dCupdate;

                }
                if (!converged)
                    break;

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
            
            // Compute the energy (interface and bulk)
            // Note: the energy is computed in the current mesh, not in the refined mesh
            auto bulk_term      = 0.25 * (1-(cnew_sol*cnew_sol).val())*(1-(cnew_sol*cnew_sol).val()); // 0.25 * (1-|c|^2)^2 double well potential
            auto grad_c         = igrad(cnew_sol, G);
            auto interface_term = 0.5  * lambda * (grad_c * grad_c.tr()).val(); // 0.5 * lambda * |grad(c)|^2

            // real_t bulk_energy       = ev.integral(meas(G) * bulk_term);
            // real_t interface_energy  = ev.integral(meas(G) * interface_term);
            // real_t total_energy      = bulk_energy + interface_energy;
            // real_t mass              = ev.integral(meas(G)*cnew_sol);

            real_t bulk_energy       = computeInt(ibasis, mp, bulk_term, A.options());
            real_t interface_energy  = computeInt(ibasis, mp, interface_term, A.options());
            real_t total_energy      = bulk_energy + interface_energy;
            real_t mass              = computeInt(ibasis, mp, cnew_sol, A.options());

            if (MESHopt.askSwitch("Adaptive",true) == 0)
            {
                csvFile << step  << ","
                        << refIt << ","
                        << A.numDofs()        << ","
                        << nSolvesStep        << ","
                        << mass               << ","
                        << bulk_energy        << ","
                        << interface_energy   << ","
                        << total_energy       << ","
                        << error_proj_c        << ","
                        << error_proj_dc       << ","
                        << assemblyTimeRefIt  << ","
                        << solverTimeRefIt    << ","
                        << projTimeRefIt << "\n";
                        csvFile.flush();  // Forces the file to write immediately
            }

            solverTimestep     += solverTimeRefIt;
            assemblyTimestep   += assemblyTimeRefIt;
            projTimestep       += projTimeRefIt;

            // Reset the errors after the loop
            error_ref_cnew  = 0;
            error_ref_cold  = 0;
            error_ref_dcnew = 0;
            error_ref_dcold = 0;

            if (MESHopt.askSwitch("Adaptive",true))
            {
                // ==================== REFINEMENT ====================
               index_t numEl = dbasis.basis(0).numElements();
                // std::vector<T> marked(numEl,false);
                ev.setIntegrationElements(dbasis);
                ev.integralElWise(meas(G) * cnew);
                std::vector<real_t> cInt = ev.elementwise();
                gsAsVector<real_t> cvec(cInt.data(),cInt.size());  // Temporary Eigen::Map
                // Compute the area of each element
                ev.integralElWise(meas(G));
                std::vector<real_t> areas = ev.elementwise();
                gsAsVector<real_t> avec(areas.data(),areas.size()); // Temporary Eigen::Map
                ev.setIntegrationElements(ibasis); //reset basis

                // Invert and normalize the element-wise average (c/area), as:
                // err = 1-|c|/a;
                cvec.array() = 1-(cvec.array().abs()/avec.array());

                clock.restart();    
                marker.setErrors(cInt);
                HElementContainer markedRef = marker.markRef();
                std::vector<index_t> refBox = marker.toRefBoxes(markedRef);
                markRefTimeRefIt += clock.stop();
                markRefTimeStep  += markRefTimeRefIt;

                size_t nElements0 =   dbasis.basis(0).numElements();
                size_t nRefined = markedRef.size();

                csvFile2 <<"," <<nRefined;
                csvFile2 <<"," <<nElements0;
                csvFile2 <<"," << static_cast<double>(nRefined)/nElements0;
                csvFile2 << '\n';                
                csvFile2.flush();        

                gsInfo<<"marked for refinement: "<<nRefined<<"\n";
                
                csvFile << step << ","
                        << refIt << ","
                        << A.numDofs() << ","
                        << nSolvesStep << ","
                        << mass << ","
                        << bulk_energy << ","
                        << interface_energy << ","
                        << total_energy << ","
                        << error_proj_c << ","
                        << error_proj_dc << ","
                        << assemblyTimeRefIt << ","
                        << solverTimeRefIt << ","
                        << projTimeRefIt << ","
                        << markRefTimeRefIt << ",";
                csvFile.flush();  // Forces the file to write immediately


                // If elements are marked for refinement
                if (nRefined!=0)
                {
                    // Refine dbasis
                    if (verbose>1) gsInfo<<"Basis before refinement:\n "<<dbasis.basis(0)<<"\n";
                    clock.restart();
                    dbasis.basis(0).refineElements(refBox);
                    ibasis.basis(0).refineElements(refBox);
                    w.setup(bc, dirichlet::l2Projection, 0);
                    refTimeRefIt += clock.stop();
                    if (verbose>1) gsInfo<<"Basis after refinement:\n "<<dbasis.basis(0)<<"\n";
                    refTimeStep       += refTimeRefIt;

                    // If there is a rule to break the refinement iteration
                    if (refIt == maxRefIt-1)
                    {
                        csvFile << refTimeRefIt << ",";
                        csvFile.flush();  // Forces the file to write immediately
                        break;  // stop refinement iteration loop if few elements are marked for refinement
                    }
                    else // next refinement iteration
                    {
                        csvFile << refTimeRefIt << "\n";
                        csvFile.flush();  // Forces the file to write immediately
                    }
                    // Add fine basis in ibasis for integration (assembler)
                    // gsInfo<<"Ibasis before: "<<ibasis.basis(0).numElements()<<"\n";
                    // ibasis.clear();
                    // for (size_t p=0; p!=dbasis.nBases(); p++)
                    //     ibasis.addBasis(dbasis.basis(p).clone());
                    // gsInfo<<"Ibasis after: "<<ibasis.basis(0).numElements()<<"\n";
                }// refine
                else
                {
                    csvFile << refTimeRefIt << ","; // to add the time of the coarsening!
                    csvFile.flush();  // Forces the file to write immediately
                    break;
                }

                // NOTE: every element weights the same, eventhough elements with larger areas can have more impact on the solution in the next refinment iteration.
                // An alternative approach would be to weight the elements by their area.
                // if (static_cast<double>(nRefined)/nElements0 < 0.02) //5%
                // {
                //     if (verbose>0) gsInfo<<"Few marked elements for refinement, stop refIt loop.\n";
                //     break; // break the refinement loop
                // }

            } // refinement switch
            else
                break; // break the refinement loop
        }// mesh adaptivity

        ibasis.clear();
        // ibasis = gsMultiBasis<>(mp_cnew); // mp_cnew has the basis from the converged RefIt
        ibasis.addBasis(dbasis.basis(0).clone()); // add the fine basis for integration

        gsInfo<< "ibasis size: "<< ibasis.totalSize() << "\n";
        gsInfo<< "dbasis size: "<< dbasis.totalSize() << "\n";


        // //! [Export visualization in ParaView]
        if (plot && step % plotmod==0)
        {
            // Export the mdesh
            collection.newTimeStep(&mp);
            collection.addField(cnew,"numerical solution");
            gsInfo << "Number of degrees of freedom:\t" << A.numDofs()  << std::endl;
            collection.saveTimeStep();
            
            if (plot_mesh)
            {
                gsMesh<> mesh(dbasis.basis(0)); // mesh refIt
                // gsMesh<> mesh_int(ibasis.basis(0)); // integration mesh
                // gsMesh<> mesh_coarse(dbasis.basis(0)); // mesh after coarsening
                gsWriteParaview(mesh, out+"/mesh_pvd/mesh_"+std::to_string(plotID),false);
                meshcollection.addPart("mesh_pvd/mesh_"+std::to_string(plotID)+".vtp",plotID,"Mesh"); //last flag is to say its a mesh
                // gsWriteParaview(mesh,"mesh_"+util::to_string(i),false);
                // meshes.addPart("mesh_"+util::to_string(i)+".vtp",i,"Mesh");
            }

            plotID++;
        }


        if (MESHopt.askSwitch("Adaptive",true))
        {

            // ==================== COARSENING ====================
            index_t numEl = dbasis.basis(0).numElements();
            ev.setIntegrationElements(dbasis);
            //     gsInfo<<"Space size = "<<w.source().size()*w.dim()<<"\n";
            //     gsInfo<<"Map size   = "<<w.mapper().mapSize()<<"\n";
            // gsInfo<< "A.numDofs() " << A.numDofs() <<"\n";
            ev.integralElWise(meas(G) * cnew);
            std::vector<real_t> cInt = ev.elementwise();
            gsAsVector<real_t> cvec(cInt.data(),cInt.size());  // Temporary Eigen::Map
            // Compute the area of each element
            ev.integralElWise(meas(G));
            std::vector<real_t> areas = ev.elementwise();
            gsAsVector<real_t> avec(areas.data(),areas.size()); // Temporary Eigen::Map
            ev.setIntegrationElements(ibasis); //reset basis

            // Invert and normalize the element-wise average (c/area), as:
            // err = 1-|c|/a;
            cvec.array() = 1-(cvec.array().abs()/avec.array());

            clock.restart();
            marker.setErrors(cInt);
            HElementContainer markedCrs = marker.markCrs(); //????????????
            std::vector<index_t> crsBox = marker.toCrsBoxes(markedCrs);
            markCrsTimeStep += clock.stop();

            size_t nCoarsened = markedCrs.size();

            // If elements are marked for refinement
            if (nCoarsened!=0)
            {
                // Refine dbasis
                if (verbose>1) gsInfo<<"Basis before coarsening:\n "<<dbasis.basis(0)<<"\n";
                clock.restart();
                // ibasis.basis(0).unrefineElements(crsBox);
                dbasis.basis(0).unrefineElements(crsBox);
                w.setup(bc, dirichlet::l2Projection, 0);
                crsTimeStep += clock.stop();
                if (verbose>1) gsInfo<<"Basis after coarsening:\n "<<dbasis.basis(0)<<"\n"; 
            }// coarsen

            csvFile << markCrsTimeStep   << ","
                    << crsTimeStep       << "\n";
            csvFile.flush();
        } // coarsening switch

        // Update time and old solutions
        solverTime     += solverTimestep;
        assemblyTime   += assemblyTimestep;
        projectionTime += projTimestep;

        // Update adaptivivity times

        markCrsTime    += markCrsTimeStep;
        crsTime        += crsTimeStep;
        markRefTime    += markRefTimeStep;
        refTime        += refTimeStep;

        time += dt_old;
        // mp_cold.swap(mp_cnew);
        // mp_dcold.swap(mp_dcnew);
        mp_cold = mp_cnew;
        mp_dcold = mp_dcnew;
    }
    if (plot)
    {
        collection.save();
        meshcollection.save();
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
    csvFile << "MarkRef  :" << markRefTime << "\n";
    csvFile << "Ref      :" << refTime << "\n";
    csvFile << "MarkCrs  :" << markCrsTime << "\n";
    csvFile << "Crs      :" << crsTime << "\n";
    csvFile << "Total     :" << assemblyTime + solverTime + projectionTime << "\n";
    csvFile << "Total (with marking):" << assemblyTime + solverTime + projectionTime + markRefTime + refTime + markCrsTime + crsTime << "\n";
    csvFile << "NSolver   :" << nSolves << "\n";
    gsInfo<<"[CLOCK] --- Time for assembly: "<<assemblyTime<<" [s]\n";
    gsInfo<<"[CLOCK] --- Time for solver  : "<<solverTime<<" [s]\n";
    gsInfo<<"[CLOCK] --- Time for projection: "<<projectionTime<<" [s]\n";
    gsInfo<<"[CLOCK] --- Number of solves: "<<nSolves<<"\n";

    csvFile.close();
    csvFile2.close();

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
    bool plot_mesh = false;
    index_t numRefine  = 1;
    index_t numElevate = 1;
    index_t verbose = 1;
    bool random = false;
    index_t projection_Crs = 0;
    index_t pattern = 0; // 0 for nucleation 1 for spinoidal decomposition
    index_t timemethod = 1; // 0 for Backward Euler, 1 for Generalized Alpha
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
    cmd.addSwitch("plotmesh", "Create a ParaView visualization file with the projection errors", plot_mesh);
    cmd.addSwitch("random", "Random initial condition of the CH problem", random);
    cmd.addInt("c", "projcoars", "Projection method for coarsening", projection_Crs);
    cmd.addInt("s", "pattern", "Phase separation pattern", pattern);
    cmd.addString( "o", "output", "Output directory", out);
    cmd.addInt("m", "timeint", "Time integration method", timemethod);

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
        solve<2,real_t>(mp, source, bc, CHopt, TIMEopt, MESHopt, Aopt, dt, maxSteps, plotmod, plot, plot_error, plot_mesh, numRefine, numElevate, verbose, random, projection_Crs, pattern, timemethod, out);
    else if (mp.geoDim()==3)
        solve<3,real_t>(mp, source, bc, CHopt, TIMEopt, MESHopt, Aopt, dt, maxSteps, plotmod, plot, plot_error, plot_mesh, numRefine, numElevate, verbose, random, projection_Crs, pattern, timemethod, out);
    else
        GISMO_ERROR("Only 2D and 3D problems are supported.");

    return EXIT_SUCCESS;

}// end main
