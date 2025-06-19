/** @file cahn-hilliard.cpp

    @brief Tutorial on how to use expression assembler to solve the Cahn-Hilliard equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Marsala (UniFi)
               H.M. Verhelst (UniFi)
               L. Venta Viñuela (UniPv)



    //  =========== QUADRATIC BASIS ===========
    ./bin/cahn-hilliard_convergence_coarsening --plot -r 6 -N 10 -c 0 -v 2 -e 1 -f pde/  .xml 

    //  =========== CUBIC BASIS ===========
    ./bin/cahn-hilliard_convergence_coarsening --plot -r 6 -N 10 -c 0 -v 2 -e 2 -f pde/  .xml 
    
    //  =========== QUARTIC BASIS ===========
    ./bin/cahn-hilliard_convergence_coarsening --plot -r 6 -N 10 -c 0 -v 2 -e 3 -f pde/  .xml 
    
    -----------------------------------------------------------------------
    TODO;
    - Change hmax to a gsExprAssembler<>::element el; el.diam();
    -----------------------------------------------------------------------

*/

//! [Include namespace]
#include <gismo.h>

using namespace gismo;
//! [Include namespace]

template <short_t dim, class T>
void solve( gsMultiPatch<T> & mp,
            gsFunctionExpr<T> & source,
            gsBoundaryConditions<T> & bc,
            gsOptionList & CHopt,
            gsOptionList & TIMEopt,
            gsOptionList & MESHopt,
            gsOptionList & Aopt,
            gsFunctionExpr<T> & ms,
            gsFunctionExpr<T> & forcingterm,
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

    out = out + "_analytical_N_" + std::to_string(maxSteps) + "_dt_" + std::to_string(dt) + "_lambda_" + std::to_string(lambda) + "_r_" + std::to_string(numRefine) + "_degree_" + std::to_string(dbasis.basis(0).maxDegree()) + "_prjCrs_" + std::to_string(projection_Crs) + "_THB_" + std::to_string(MESHopt.askSwitch("THB")) + "_Adaptive_" + std::to_string(MESHopt.askSwitch("Adaptive")) + "_refIt_" + std::to_string(MESHopt.askInt("RefIt"));

    //! [Prepare the basis]

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    A.options().setSwitch("SameElement",Aopt.askSwitch("SameElement",true));
#ifdef GISMO_WITH_PARDISO
    A.options().addString("LinearSolver", "Name of the linear solver to be used", "PardisoLU");
#   else
    A.options().addString("LinearSolver", "Name of the linear solver to be used", "SimplicialLDLT");
#endif

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
    ev.options().update(A.options(),gsOptionList::addIfUnknown); 
    ev.setIntegrationElements(dbasis);

    gsDebugVar(dbasis.basis(0).maxDegree());
    gsDebugVar(ibasis.basis(0));

    // Set the geometry map
    geometryMap G = A.getMap(mp);

    // Set the discretization space
    space w = A.getSpace(dbasis);

    // basis.init(dbasis, cf);

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

    gsMultiPatch<> mp_cold, mp_cnew;
    auto     cold  = A.getCoeff(mp_cold );
    auto     cnew  = A.getCoeff(mp_cnew );
    // %%%%%%%%%%%%%%%%%%%%%%%% Analytical intial condition %%%%%%%%%%%%%%%%%%%%%%%%
    GISMO_ASSERT(mp.geoDim()==ms.domainDim(),"Domain dimension of the source function should be equal to the geometry dimension, but "<<ms.domainDim()<<"!="<<mp.geoDim());
    gsMatrix<> tmp;
    gsInfo << ms << "\n";
    ms.set_t(0); // set time to zero in manufactured solution (initial condition)
    real_t error = gsL2Projection<real_t>::project(dbasis,mp,ms,tmp,A.options());  // 3rd arg has to be multipatch
    if (verbose>0) gsInfo << "L2 projection error "<<error<<"\n";
    mp_cold.addPatch(dbasis.basis(0).makeGeometry(tmp));
    mp_cnew = mp_cold;

    real_t Q0norm = 1, Qnorm = 10;
    real_t tol = TIMEopt.askReal("tol",1e-4);

    gsParaviewCollection collection(out+"/solution", &ev);
    collection.options().setSwitch("plotElements", true);
    collection.options().setInt("plotElements.resolution", 4);
    collection.options().setInt("numPoints",(mp.geoDim()==3) ? 10000 : 5000);
    collection.options().setInt("precision", 10); // digits 10^-10

    real_t time = 0;
    bool converged = false;

    // Sparse matrix for Nitsche contribution
    gsSparseMatrix<> K_nitsche; // empty variable

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
        ms.set_t(0); // set time to zero in manufactured solution (initial condition)
        auto u_manufactured = ev.getVariable(ms, G);
        collection.addField(u_manufactured,"analytical solution");
        collection.saveTimeStep();
    }

    if (plot)
        collection.save(); 

    // =================== Save the coarsened basis !!! ===================
    gsFileData<> fdbas;
    fdbas.add(dbasis,0);
    fdbas.save("basis_cir_"+std::to_string(MESHopt.getReal("CoarsenParam"))+"_r_"+std::to_string(numRefine)+"_d_"+std::to_string(dbasis.basis(0).maxDegree())+".xml");
    gsInfo << "Exported to "<<  "basis_cir_"+std::to_string(MESHopt.getReal("CoarsenParam"))+"_r_"+std::to_string(numRefine)+"_d_"+std::to_string(dbasis.basis(0).maxDegree())+".xml" <<"\n";
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

    gsFunctionExpr<> ms;
    fd.getId(7, ms); // id=7: manufactured solution

    gsFunctionExpr<> forcingterm;
    fd.getId(8, forcingterm); // id=8: forcing term

    if (mp.geoDim()==2)
        solve<2,real_t>(mp, source, bc, CHopt, TIMEopt, MESHopt, Aopt, ms, forcingterm, dt, maxSteps, plotmod, plot, plot_error, numRefine, numElevate, verbose, random, projection_Crs, pattern, out);
    else if (mp.geoDim()==3)
        solve<3,real_t>(mp, source, bc, CHopt, TIMEopt, MESHopt, Aopt, ms, forcingterm, dt, maxSteps, plotmod, plot, plot_error, numRefine, numElevate, verbose, random, projection_Crs, pattern, out);
    else
        GISMO_ERROR("Only 2D and 3D problems are supported.");

    return EXIT_SUCCESS;

}// end main
