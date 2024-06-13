/** @file cahn-hilliard.cpp

    @brief Tutorial on how to use expression assembler to solve the Cahn-Hilliard equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Marsala (UniFi)
               H.M. Verhelst (UniFi)
               L. Venta Viñuela (UniPv)
    
    
    Run a simple Cahn-Hilliard example with an analytical initial condition "0.1 * cos(2*pi*x) * cos(2*pi*y)" (Strong enforcement) (Gomez et al., 2014)
    ./bin/cahn-hilliard_example --plot -N 80 --plot
    
    Run a simple Cahn-Hilliard example with an analytical initial condition "0.1 * cos(2*pi*x) * cos(2*pi*y)" (Nitsche) (Bracco et al., 2023)
    ./bin/cahn-hilliard_example --plot -N 80 --nitsche --plot

    Run a simple Cahn-Hilliard example with a random normal initial concentration distribution of mean 0.0 until (almost) equilibrium (Nitsche)
    ./bin/cahn-hilliard_example --plot -N 1000 --nitsche --initial --plot
    

    -----------------------------------------------------------------------
    TODO;
    - Change hmax to a gsExprAssembler<>::element el; el.diam();
    -----------------------------------------------------------------------



*/

//! [Include namespace]
#include <gismo.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    real_t dt = 1e-3;
    index_t maxSteps = 10;

    index_t plotmod = 1;

    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 1;
    index_t numElevate = 1;

    index_t verbose = 1;
    bool random = false;

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
    cmd.addSwitch("random", "Random initial condition of the CH problem", random);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]
    gsFileData<> fd(fn);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    gsMultiPatch<> mp;
    fd.getId(0, mp); // id=0: Multipatch domain

    gsFunctionExpr<> source;
    fd.getId(1, source); // id=1: initial condition function
    gsInfo<<"Initial condition function "<< source << "\n";

    gsBoundaryConditions<> bc;
    fd.getId(2, bc); // id=2: boundary conditions
    bc.setGeoMap(mp);
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";

    gsOptionList CHopt;
    fd.getId(3, CHopt); // id=3: reference solution

    real_t theta    = CHopt.askReal("theta",1.5);
    real_t lambda   = CHopt.askReal("lambda",1/(32*pow(EIGEN_PI,2)));
    real_t M0       = CHopt.askReal("M0",0.005);
    real_t penalty  = CHopt.askReal("penalty",1e4*lambda);

    gsOptionList TIMEopt;
    fd.getId(4, TIMEopt); // id=4: time integrator options

    gsOptionList Aopt;
    fd.getId(5, Aopt); // id=5: assembler options

    gsOptionList MESHopt;
    fd.getId(6, MESHopt); // id=6: mesher options
    //! [Read input file]

    //! [Prepare the basis]
    gsMultiBasis<> dbasis_tmp(mp,true);
    gsMultiBasis<> dbasis;

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis_tmp.setDegree( dbasis_tmp.maxCwiseDegree() + numElevate);


    gsDebugVar(dbasis_tmp);
    // Cast every basis of dbasis to a gsTHBSplineBasis
    for (size_t p=0; p!=dbasis_tmp.nBases(); p++)
    {
        // TODO: Make dimension-independent over the template
        if (gsTensorBSplineBasis<2,real_t> * b = dynamic_cast<gsTensorBSplineBasis<2,real_t>*>(&dbasis_tmp.basis(p)))
            dbasis.addBasis(new gsTHBSplineBasis<2,real_t>(*b));
        else if (gsTHBSplineBasis<2,real_t> * b = dynamic_cast<gsTHBSplineBasis<2,real_t>*>(&dbasis_tmp.basis(p)))
            dbasis.addBasis(b->clone());
        else
            GISMO_ERROR("Basis is neither a gsTHBSplineBasis nor a gsTensorBSplineBasis");

        // Refine the basis for `numRefine` levels
        gsMatrix<> box = dbasis.basis(p).support();
        for (index_t r = 0; r!=numRefine; r++)
            dbasis.basis(p).refine(box);
    }




    // Determine maximum mesh size
    real_t hmax = 0;
    for (size_t p=0; p!=dbasis.nBases(); p++)
        hmax = math::max(hmax, dbasis.basis(p).getMaxCellLength());

    for (size_t p=0; p!=dbasis.nBases(); p++)
        gsDebug<<dbasis.basis(p).maxDegree()<<"\n";

    //! [Prepare the basis]

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //! [Problem setup]
    gsExprAssembler<> A(1,1);

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);

    // Set the geometry map
    // geometryMap G = A.getMap(surface);
    geometryMap G = A.getMap(mp);

    // Set the discretization space
    space w = A.getSpace(dbasis);

    // basis.init(dbasis, cf);

    // Solution vector and solution variable
    gsMatrix<> Cnew, Calpha, Cold;
    gsMatrix<> dCnew,dCalpha,dCold, dCupdate;

    // Solution variables for the intermediate solutions (during time integration)
    solution c = A.getSolution(w, Calpha); // C
    solution dc = A.getSolution(w, dCalpha); // \dot{C}

    // Solution variables for the previous and next solutions (before and after time step)
    solution cold = A.getSolution(w, Cold); // Cold
    solution cnew = A.getSolution(w, Cnew); // Cnew
    solution dcold = A.getSolution(w, dCold); // \dot{Cold}
    solution dcnew = A.getSolution(w, dCnew); // \dot{Cnew}

    // Derivatives of the double well potential (Gomez et al., 2008)
    auto dmu_c = - 1.0 + 3.0 * (c*c).val(); // f_2 (second derivative of double well)
    auto ddmu_c = 6*c.val(); // f_3 (third derivative of double well)

    // Mobility
    auto M_c  = 1.0 + 0.0*c.val(); // replace with const_expr(1.0) instead of using 0*c
    auto dM_c = 0.0 * igrad(c,G); // replace with const_expr(1.0) instead of using 0*c!!

    auto residual = w*dc + // M
                    M_c.val() * igrad(w,G)  * dmu_c * igrad(c,G).tr() + // F_bar
                    M_c.val() * ilapl(w,G)*lambda*ilapl(c,G).val(); // K_laplacian
                    // lambda*ilapl(c,G).val()*igrad(w,G)*dM_c.tr() + // term gradient mobility!

    //! [Problem setup]

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

    // Generalized-alpha method parameters
    real_t rho_inf = TIMEopt.askReal("rho_inf",0.5);
    real_t alpha_m = 0.5*(3-rho_inf) / (1+rho_inf);
    real_t alpha_f = 1 / (1+rho_inf);
    real_t gamma   = 0.5 + alpha_m - alpha_f;
    // time stepping options
    index_t maxIt = 50;

    gsMatrix<> Q, Q1, Q2;
    gsSparseMatrix<> K, K_m, K_f;

    // Legend:
    // C_old   = C_n
    // C       = C_n+1,i-1
    // C_alpha = C_{n+alpha_f,i-1}
    // dC

    gsInfo<<"Starting.."<<"\n";



    gsInfo<<"Initial condition.."<<"\n";

    if (random)
    {
        // %%%%%%%%%%%%%%%%%%%%%%%% Random initial condition %%%%%%%%%%%%%%%%%%%%%%%%
        gsMatrix<> tmp = gsMatrix<>::Random(A.numDofs(),1);
        Cold = tmp.array()*CHopt.askReal("ampl",0.005); //random uniform variable in [-0.05,0.05]
        Cold.array() += CHopt.askReal("mean",0.0); // 0.45
    }
    else
    {
        // %%%%%%%%%%%%%%%%%%%%%%%% Analytical intial condition %%%%%%%%%%%%%%%%%%%%%%%%
        GISMO_ASSERT(mp.geoDim()==source.domainDim(),"Domain dimension of the source function should be equal to the geometry dimension, but "<<source.domainDim()<<"!="<<mp.geoDim());
        gsMatrix<> tmp;
        Cold.setZero(A.numDofs(),1);
        real_t error = gsL2Projection<real_t>::projectFunction(dbasis,source,mp,tmp);  // 3rd arg has to be multipatch
        if (verbose>0) gsInfo << "L2 projection error "<<error<<"\n";
        for (index_t i = 0; i < dbasis.basis(0).size(); i++)
            if (w.mapper().is_free(i))
                Cold(w.mapper().index(i),0) = tmp(i,0);
    }

    Calpha = Cold;
    dCold.setZero(A.numDofs(),1);

    real_t Q0norm = 1, Qnorm = 10;
    real_t tol = TIMEopt.askReal("tol",1e-4);

    gsParaviewCollection collection("ParaviewOutput/solution", &ev);
    collection.options().setSwitch("plotElements", true);
    collection.options().setInt("plotElements.resolution", 4);
    collection.options().setInt("numPoints",(mp.geoDim()==3) ? 10000 : 1000);

    real_t dt_old = dt;
    real_t t_rho = TIMEopt.askReal("t_rho",0.9);
    real_t t_err = 1;
    index_t lmax = 1;
    std::vector<gsMatrix<>> Csols(2);

    real_t tmp_alpha_m = 1;
    real_t tmp_alpha_f = 1;
    real_t tmp_gamma   = 1;

    real_t time = 0;
    bool converged = false;

    // Sparse matrix for Nitsche contribution
    gsSparseMatrix<> K_nitsche; // empty variable


    // ! [Load mesher options]
    // DIMENSION-INDEPENDENT
    // gsAdaptiveMeshingBase<real_t> * mesher;
    // if(mp.geoDim()==2)
    //     mesher = new gsAdaptiveMeshing<2,real_t>(dbasis);
    // else if(mp.geoDim()==2)
    //     mesher = new gsAdaptiveMeshing<3,real_t>(dbasis);
    gsAdaptiveMeshing<2,real_t> mesher(dbasis);
    mesher.options().setInt("RefineRule",MESHopt.askInt("RefineRule",1));
    mesher.options().setInt("CoarsenRule",MESHopt.askInt("CoarsenRule",1));
    mesher.options().setReal("RefineParam",MESHopt.askReal("RefineParam",0.1));
    mesher.options().setReal("CoarsenParam",MESHopt.askReal("CoarsenParam",0.1));
    mesher.options().setSwitch("Admissible",MESHopt.askSwitch("Admissible",false));
    mesher.options().setInt("MaxLevel",numRefine);
    mesher.options().setSwitch("Absolute",MESHopt.askSwitch("Absolute",true));
    mesher.getOptions();
    // ! [Load mesher options]



    for (index_t step = 0; step!=maxSteps; step++)
    {
        for (index_t refIt = 0; refIt!=MESHopt.askInt("RefIt",5); refIt++)
        {
            // Resize the data structure inside the mesher
            mesher.rebuild();

            // Reset the assembler
            A.initSystem();

            // Setup the space (compute Dirichlet BCs)
            w.setup(bc, dirichlet::l2Projection, 0);

            for (index_t dt_it = 0; dt_it != lmax; dt_it++)
            {
                if (verbose>0) gsInfo<<"Time step "<<step<<"/"<<maxSteps<<", iteration "<<dt_it<<": dt = "<<dt<<", [t_start,t_end] = ["<<time<<" , "<<time+dt<<"]"<<"\n";
                tmp_alpha_m = tmp_alpha_f = tmp_gamma = 1;

                for (index_t k = 0; k!=2; k++)
                {
                    converged = false;
                    std::string method = (k==0) ? "Backward Euler " : "Generalized Alpha ";
                    // ==================================================================================
                    // Predictor
                    Cnew = Cold;
                    dCnew = (tmp_gamma-1)/tmp_gamma * dCold;

                    Q0norm = 1;
                    Qnorm = 10;

                    for (index_t it = 0; it!= maxIt; it++)
                    {
                        A.initMatrix();
                        A.clearRhs(); // Resets to zero the values of the already allocated to residual (RHS)
                        Calpha.noalias()  = Cold  + tmp_alpha_f * ( Cnew  - Cold );
                        dCalpha.noalias() = dCold + tmp_alpha_m * ( dCnew - dCold);

                        // Assemble the RHS
                        A.assemble(residual * meas(G));
                        Q = A.rhs();

                        // Assemble the Nitsche BC on the sides with Neumann condition
                        // A.clearMatrix(); // Resets to zero the values of the already allocated to matrix (LHS)
                        A.initMatrix();
                        A.assembleBdr(bc.get("Neumann"), - lambda * igrad(w,G) *  nv(G)  * ilapl(w,G).tr() + // consistency term
                                      penalty * (igrad(w,G) * nv(G).normalized()) * hmax * (igrad(w,G) * nv(G)).tr() - // penalty (stabilizing) term
                                      lambda * ilapl(w,G) * (igrad(w,G)  * nv(G)).tr()); // symmetry term
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
                        A.assemble(meas(G) * (w*w.tr()*tmp_alpha_m +// K_m
                                            (tmp_alpha_f * tmp_gamma * dt)* (dmu_c *igrad(w,G) * igrad(w,G).tr() + // K_f1
                                            ddmu_c * igrad(w,G) * igrad(c,G).tr() * w.tr() + // K_f2
                                            lambda * ilapl(w,G) * ilapl(w,G).tr()))); // K_laplacian
                                            // lambda * igrad(w,G)*dM_c.tr()*ilapl(w,G).tr()   +  // K_mobility

                        K = A.giveMatrix();
                        if (bc.get("Neumann").size()!=0)
                            K += (tmp_alpha_f * tmp_gamma * dt) * K_nitsche; // add the Nitsche term to the stiffness matrix

    #ifdef GISMO_WITH_SUPERLU
                        if (0==k)
                            solver.analyzePattern(K);
                        solver.factorize(K);
    #else
                        solver.compute(K);
    #endif
                        dCupdate = solver.solve(-Q);

                        dCnew += dCupdate;
                        Cnew.noalias() += (tmp_gamma*dt)*dCupdate;
                    }
                    if (!converged)
                        break;

                    // %% Switch to generalized-alpha parameters (k=1)
                    tmp_alpha_m = alpha_m;
                    tmp_alpha_f = alpha_f;
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
            // -------------REFINEMENT-------------------
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

            // Mark the elements for refinement
            gsHBoxContainer<2,real_t> refine;//, coarsen;
            mesher.markRef_into(cInt,refine);

            // If elements are marked for refinement
            if (refine.totalSize()!=0)
            {
                // Store the old and new solutions on the previous (SINGLE PATCH ASSUMPTION)
                // Take their full coefficient vector (including eliminated DoFs)
                gsMatrix<> CnewF, dCnewF, ColdF, dColdF;
                cold.extractFull(ColdF);
                dcold.extractFull(dColdF);
                cnew.extractFull(CnewF);
                dcnew.extractFull(dCnewF);

                // Create the geometry objects
                gsGeometry<>::uPtr Cold_ = dbasis.basis(0).makeGeometry(give(ColdF));
                gsGeometry<>::uPtr dCold_ = dbasis.basis(0).makeGeometry(give(dColdF));
                gsGeometry<>::uPtr Cnew_ = dbasis.basis(0).makeGeometry(give(CnewF));
                gsGeometry<>::uPtr dCnew_ = dbasis.basis(0).makeGeometry(give(dCnewF));

                // Refine dbasis
                if (verbose>1) gsInfo<<"Basis before refinement:\n "<<dbasis.basis(0)<<"\n";
                mesher.refine(refine);
                if (verbose>1) gsInfo<<"Basis after refinement:\n "<<dbasis.basis(0)<<"\n";

                // Project the old and new solutions onto the new basis
                gsQuasiInterpolate<real_t>::localIntpl(dbasis.basis(0),*Cold_,ColdF);
                gsQuasiInterpolate<real_t>::localIntpl(dbasis.basis(0),*dCold_,dColdF);
                gsQuasiInterpolate<real_t>::localIntpl(dbasis.basis(0),*Cnew_,CnewF);
                gsQuasiInterpolate<real_t>::localIntpl(dbasis.basis(0),*dCnew_,dCnewF);

                // Setup the space to obtain a new DoF mapper
                w.setup(bc, dirichlet::l2Projection, 0);
                // Resize the new solution vectors (which exclude eliminated DoFs)
                Cold.resize(w.mapper().freeSize(),1);
                dCold.resize(w.mapper().freeSize(),1);
                Cnew.resize(w.mapper().freeSize(),1);
                dCnew.resize(w.mapper().freeSize(),1);

                // Insert the interpolated coefficients inside the solution objects
                for (index_t i = 0; i < dbasis.basis(0).size(); i++)
                    if (w.mapper().is_free(i))
                    {
                        Cold(w.mapper().index(i),0) = ColdF(i,0);
                        dCold(w.mapper().index(i),0) = dColdF(i,0);
                        Cnew(w.mapper().index(i),0) = CnewF(i,0);
                        dCnew(w.mapper().index(i),0) = dCnewF(i,0);
                    }
            }// refine
            else
                break;

        }// mesh adaptivity

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
        gsHBoxContainer<2,real_t> coarsen;
        mesher.markCrs_into(cInt,coarsen);

        // If elements are marked for refinement
        if (coarsen.totalSize()!=0)
        {
            // Store the old and new solutions on the previous (SINGLE PATCH ASSUMPTION)
            // Take their full coefficient vector (including eliminated DoFs)
            gsMatrix<> CnewF, dCnewF, ColdF, dColdF;
            cnew.extractFull(CnewF);
            dcnew.extractFull(dCnewF);

            // Create the geometry objects
            gsGeometry<>::uPtr Cnew_ = dbasis.basis(0).makeGeometry(give(CnewF));
            gsGeometry<>::uPtr dCnew_ = dbasis.basis(0).makeGeometry(give(dCnewF));

            // Refine dbasis
            if (verbose>1) gsInfo<<"Basis before coarsening:\n "<<dbasis.basis(0)<<"\n";
            mesher.unrefine(coarsen);
            if (verbose>1) gsInfo<<"Basis after coarsening:\n "<<dbasis.basis(0)<<"\n";

            // Project the old and new solutions onto the new basis
            gsQuasiInterpolate<real_t>::localIntpl(dbasis.basis(0),*Cnew_,CnewF);
            gsQuasiInterpolate<real_t>::localIntpl(dbasis.basis(0),*dCnew_,dCnewF);

            // Setup the space to obtain a new DoF mapper
            w.setup(bc, dirichlet::l2Projection, 0);
            // Resize the new solution vectors (which exclude eliminated DoFs)
            Cnew.resize(w.mapper().freeSize(),1);
            dCnew.resize(w.mapper().freeSize(),1);

            // Insert the interpolated coefficients inside the solution objects
            for (index_t i = 0; i < dbasis.basis(0).size(); i++)
                if (w.mapper().is_free(i))
                {
                    Cnew(w.mapper().index(i),0) = CnewF(i,0);
                    dCnew(w.mapper().index(i),0) = dCnewF(i,0);
                }
        }// coarsen


        // Update time and old solutions
        time += dt_old;
        Cold = Cnew;
        dCold = dCnew;

        //! [Export visualization in ParaView]
        if (plot && step % plotmod==0)
        {
            // Export the mesh
            collection.newTimeStep(&mp);
            collection.addField(cnew,"numerical solution");

            collection.saveTimeStep();
        }
    }

    if (plot)
    {
        collection.save();
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";

    return EXIT_SUCCESS;

}// end main
