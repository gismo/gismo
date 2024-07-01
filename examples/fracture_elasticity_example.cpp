/** @file cahn-hilliard.cpp

    @brief Tutorial on how to use expression assembler to solve the Cahn-Hilliard equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): 
    
*/

//! [Include namespace]
#include <gismo.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    real_t theta = 1.5;
    real_t lambda = 1/(32*pow(EIGEN_PI,2));
    // real_t lambda = 6.15e-4;
    real_t L0 = 1;
    real_t M0 = 0.005;
    real_t dt = 1e-3;
    index_t maxSteps = 10;

    // real_t eps_penalty = 1e4 *lambda;   // not sure about the value!


    bool output = true;


    index_t plotmod = 1;

    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 1;
    index_t numElevate = 1;
    bool last = false;
    bool nitsche = false;
    real_t mean = 0.0;

    bool random = false;

    bool do3D = false;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addReal( "t", "dt","dt parameter",dt); // -t () or --dt ()
    cmd.addReal( "T", "theta","Theta parameter",theta);
    cmd.addReal( "l", "lambda","lambda parameter",lambda);
    cmd.addReal( "L", "L0","L0 parameter",L0);
    cmd.addReal( "M", "M0","M0 parameter",M0);
    cmd.addInt ( "N", "Nsteps", "Number of time steps",  maxSteps );
    cmd.addInt ( "p", "PlotMod", "Modulo for plotting",  plotmod );
    // cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("nitsche", "Weak BC enforcement with Nitsche", nitsche);
    cmd.addReal("m","mean", "Mean value of the normal random initial condition", mean);
    cmd.addSwitch("initial", "Initial random condition of the CH problem", random);
    cmd.addSwitch("3D", "Flag for 3D CH", do3D);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    // %%%%%%%%%%%%%%%%%% Definition of the geometry and the basis %%%%%%%%%%%%%%%%%%
    // 0. Single patch construction parameters
    index_t n = 20;
    index_t degree = 2;
    real_t hmax = 1.0/(n-degree);

    // 1. construction of a knot vector for each direction
    // n - degree - 1 interior knots
    // degree + 1 multiplicity at the ends
    // n-p elements!
    gsKnotVector<> kv(0, 1, n - degree - 1, degree + 1); // check definition of knot vectors

    // 2. construction of a basis    
    gsMultiBasis<> dbasis;
    if (do3D)
        dbasis.addBasis(new gsTensorBSplineBasis<3>(kv, kv, kv) );
    else
        dbasis.addBasis(new gsTensorBSplineBasis<2>(kv, kv) );

    // 3. Construction of a square or cube
    gsMultiPatch<> mp;
    if (do3D)
        mp.addPatch( *gsNurbsCreator<>::BSplineCube(1) );
    else
        mp.addPatch( *gsNurbsCreator<>::BSplineSquare(1) );
    mp.computeTopology();
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    // Boundary conditions
    gsBoundaryConditions<> bc;
    // ============================== Add Dirichlet BC ==============================
    gsVector<> displ_bc (2); //could be a scalar? check!
    displ_bc<< 0,1; // 0 in x 1(quantity) in y
    gsConstantFunction<> displ_bc_fun (displ_bc,2);
    bc.addCondition(boundary::side::north,condition_type::dirichlet,&displ_bc_fun,0,false,1); //last componnent is the component x(0), y(1), z(2)
    bc.addCondition(boundary::side::south,condition_type::dirichlet,0,0,false,-1); //last componnent -1 fixes all displacements
    // BC for phase field (fracture)
    // bc.addCondition(boundary::side::north,condition_type::dirichlet,&displ_bc_fun,1); // 1 because of the ID of the space (phase field is 1)
    bc.setGeoMap(mp);    
    // ==============================================================================


    //! [Problem setup]
    gsExprAssembler<> A(2,2); //2 trial 2 test spaces

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);

    // Set the geometry map
    geometryMap G = A.getMap(mp);

    // Set the discretization space (!!!!!!!!!!different for displacements and phase field!!!!!!!!!)
    space w = A.getSpace(dbasis,2,0); // displacements (3rd argument is the id of the space)
    space b = A.getSpace(dbasis,1,1); // phase field d

    // Solution vectors and solution variables
    gsMatrix<> Unew, Uold;
    gsMatrix<> Dnew, Dold;
    
    solution uold = A.getSolution(w, Uold); 
    solution unew = A.getSolution(w, Unew); 

    solution dold = A.getSolution(w, Dold); 
    solution dnew = A.getSolution(b, Dnew); 

    //! [Problem setup]

    // Define linear solver (install SuperLUMT-devel)
#ifdef GISMO_WITH_SUPERLU
    gsSparseSolver<>::SuperLU solver;
#   else
    gsSparseSolver<>::LU solver;
#endif
    
    // Generalized-alpha method parameters
    // real_t rho_inf = 0.5;
    // real_t alpha_m = 0.5*(3-rho_inf) / (1+rho_inf);
    // real_t alpha_f = 1 / (1+rho_inf);
    // real_t gamma   = 0.5 + alpha_m - alpha_f;

    // Legend:
    // C_old   = C_n
    // C       = C_n+1,i-1
    // C_alpha = C_{n+alpha_f,i-1}
    // dC

    gsInfo<<"Starting.."<<"\n";
    
    // Setup the spaces (compute Dirichlet BCs)
    w.setup(bc, dirichlet::l2Projection, 0); // computes the values on CPs of dirichlet BC   
    b.setup(bc, dirichlet::l2Projection, 0);

    gsInfo<<"Initial condition.."<<"\n";

    // AT2 parameters                                        
    real_t c_w = 2.0;
    // Material properties
    real_t G_c = 1e4; // in N/m
    real_t ell = 0.015e-3; // in m
    real_t E = 1e11; // Young's modulus
    real_t nu = 0.25; // Poisson's ratio
    // Formulas 3D
    real_t mu = E/(2.0*(1+nu)); // shear modulus 
    real_t kappa = E/(3.0*(1-2.0*nu)); // bulk modulus

    // ================ Initialize displacements and phase field ================
    Dold.setZero(A.numDofs(),1);
    Uold.setZero(A.numDofs(),1);

    // 1) Knot coordinates x and y ?? if not in parametric space?
    // 2) Check if values of x and y are within the bounds [0, 7.5e-3] and [0.5-l0/3, 0.5+l0/3], respectively
    // 3) Assign value of 1 
    // 4) L2 project to the control points (how do I L2 project without a source function)
    // gsL2Projection<real_t>::projectFunction(dbasis,source,mp,tmp); 

    gsFunctionExpr<> fun("if (( 0 < x < 7.5e-3) and (0.5-0.005e-3 < y < 0.5+0.005e-3)) { 1 } else { 0 };",2);
    gsGeometry<>::uPtr geom_ptr = dbasis.basis(0).makeGeometry(give(coefs(dbasis.size(),2)));        
    gsWriteParaview(*geom_ptr,mp,"fun");
    gsMatrix<> tmp;
    gsL2Projection<real_t>::projectFunction(dbasis,fun,mp,tmp); 
    //if (verbose>0) gsInfo << "L2 projection error "<<error<<"\n";
    for (index_t i = 0; i < dbasis.basis(0).size(); i++)
    if (b.mapper().is_free(i))
        Dold(b.mapper().index(i),0) = tmp(i,0);
    
    Dnew = Dold;
    Unew = Uold;
    // ============================================================================

    // Loop variables
    real_t res_u_norm0, res_u_norm, res_d_norm0, res_d_norm;
    // real_t res_norm0 = 1, res_norm = 10;
    real_t tol_NR = 1e-6;
    real_t tol_stag = 1e-4;
    index_t maxStag = 10; // Staggered iteration
    index_t maxIt = 50; // NR iteration
    real_t penalty_irrev = 2.7e12;


    gsParaviewCollection collection("ParaviewOutput/solution", &ev);
    collection.options().setSwitch("plotElements", true);
    collection.options().setInt("plotElements.resolution", 4);
    collection.options().setInt("numPoints",(do3D) ? 10000 : 1000);

    // real_t dt_old = dt;
    // real_t t_rho = 0.9;
    // real_t t_err = 1;
    // index_t lmax = 1;
    // std::vector<gsMatrix<>> Csols(2);

    // real_t tmp_alpha_m = 1;
    // real_t tmp_alpha_f = 1;
    // real_t tmp_gamma   = 1;

    real_t time = 0;
    bool converged = false;

    A.initSystem(); // Initialize the system (outside the loops)
    
    // In 2D (?)
    // gsMatrix<> P_V(3,3);
    // gsMatrix<> P_D(3,3);

    gsSparseMatrix<> K_u, K_d;
    gsMatrix<> Q_u, Q_d;

    // Trace and deviatoric matrices (?)
    // P_V << 1,1,0,  1,1,0,  0,0,0;
    // P_D << 2.0/3.0,-1.0/3.0,0,  -1.0/3.0,2.0/3.0,0,  0,0,1.0/2.0;

    gsMatrix<> P_V(9,1);
    gsMatrix<> P_D(9,1);

    P_V.reshape(3,3)<< 1,1,0,  1,1,0,  0,0,0;
    P_D.reshape(3,3)<< 2.0/3.0,-1.0/3.0,0,  -1.0/3.0,2.0/3.0,0,  0,0,1.0/2.0;

    gsConstantFunction<> PVmat(P_V,2);
    gsConstantFunction<> PDmat(P_D,2);

    auto PV = A.getCoeff(PVmat);
    auto PD = A.getCoeff(PDmat);

 
    // I need to do it inside the loop (?)
    auto eps = 0.5*(igrad(unew,G).tr() + igrad(unew,G)); // linearized strain  //  I want the strains! (size 2x2 or 3X1 in 2D (plane strain))
    auto eps_v = eps.trace(); //volumetric strain (scalar) 
    
    // ============ Needs correction bc not an expression but a variable! ============
    auto Heaviside_pos = ternary(eps_v, 1, 0); // positive Heaviside function
    auto Heaviside_neg = ternary(-eps_v, 1, 0); // negative Heaviside function
    // auto Heaviside_pos = (eps_v > 0) ? 1 : 0;
    // auto Heaviside_neg = (-eps_v > 0) ? 1 : 0; 
    // ===============================================================================
    // auto D_pos = kappa * Heaviside_pos * P_V + 2*mu*P_D;
    // auto D_neg = kappa * Heaviside_neg * P_V; // 3x3 matrix

    auto g_d = pow(1-dnew.val(),2) + 0.0; // degradation function (Add term to prevent numerical difficulties when d approaches 1)

    // umax =  0.04 mm
    // delta_u = 5e-3 mm


    for (index_t step = 0; step!=maxSteps; step++)
    {
        // gsInfo<<"Time step "<<step<<"/"<<maxSteps<<", iteration "<<dt_it<<": dt = "<<dt<<", [t_start,t_end] = ["<<time<<" , "<<time+dt<<"]"<<"\n";
        gsInfo<<"Time step "<<step<<"/"<<maxSteps<<"\n";
        res_u_norm0 = 1;
        res_d_norm0 = 1;
        res_u_norm = 10;
        res_d_norm = 10;

        // Increase the applied displacement!

        for (index_t stag = 0; stag!=maxStag; stag++)
        {
            for (index_t it_eq = 0; it_eq!=maxIt; it_eq++) 
            {
                // Pending: Dirichlet boundary condition!
                A.clearMatrix(); //?
            
                auto D_pos = kappa * Heaviside_pos * reshape(PV,3,3) + 2*mu*reshape(PD,3,3);
                auto D_neg = kappa * Heaviside_neg * reshape(PV,3,3); // 3x3 matrix

                auto D_split_dmg = g_d * D_pos + D_neg;// material matrix with energy split (linear elasticity)

                K_u = igrad(w,G) * D_split_dmg * igrad(w,G).tr(); // stiffness matrix elastic problem

                A.assemble(K_u * meas(G)); // stiffness matrix
                K_u = A.matrix();
                solver.compute(K_u);

                // Boundary term (Body forces? add missing terms)
                auto f = 0; //change this
                A.assembleBdr(bc.get("Neumann"), w * f * nv(G)); // missing the f function!
                Q_u = A.rhs();
                
                // Compute solution
                Unew = solver.solve(Q_u);
                
                // Compute residual norm
                real_t res_u =  (K_u * Unew - Q_u).norm();

                // Convergence check
                if (it_eq == 0) res_u_norm0 = res_u;
                else         res_u_norm = res_u;
                gsInfo<<"\t\tNR EQ iter   "<<it_eq<<": res_u = "<<res_u_norm/res_u_norm0<<"\n";

                if (it_eq>0 && res_u_norm/res_u_norm0 < tol_NR)
                {
                    gsInfo<<"\t\tEquilibrium converged in "<<it_eq<<" iterations\n";
                    converged = true;
                    break;
                }
                else if (it_eq==maxIt-1)
                {
                    gsInfo<<"\t\tEquilibrium did not converge!\n";
                    converged = false;
                    break;
                }

            } // NR equilibrium equation

            for (index_t it_pf = 0; it_pf!=maxIt; it_pf++)
            {
                auto strain_en_pos = 0.5 * igrad(w,G) * D_pos * igrad(u,G); // (scalar!!! check dimensions) positive part of the strain energy
                
                A.clearMatrix(); 
                
                real_t penalty = 0.0;
                // if alpha - alpha_old is negative, change penalty parameter to take into account integral in the formulation!
                if (dnew < dold)
                    penalty = penalty_irrev;

                // if (Dnew - Dold < 0)
                //     penalty = penalty_irrev;

                
                // Stiffness matrix
                K_d = (w*w.tr()) * (2.0 * strain_en_pos + G_c/c_w * (2.0/ell) + penalty)  +   
                      (G_c * 2.0 * ell / c_w) * igrad(w,G)* igrad(w,G).tr();
                A.assemble(K_d * meas(G)); // stiffness matrix
                K_d = A.matrix();
                solver.compute(K_d);

                // RHS
                A.assemble(2.0 * w * strain_en_pos + penalty * w*w.tr() * meas(G));
                Q_d = A.rhs();

                // Compute solution
                Dnew = solver.solve(Q_d);
                
                // Compute residual norm
                real_t res_d =  (K_d * Dnew - Q_d).norm();
                
                // Convergence check
                if (it_pf == 0) res_d_norm0 = res_d;
                else         res_d_norm = res_d;
                gsInfo<<"\t\tNR PF iter   "<<it_pf<<": res_d = "<<res_d_norm/res_d_norm0<<"\n";

                if (it_pf>0 && res_d_norm/res_d_norm0 < tol_NR)
                {
                    gsInfo<<"\t\tPF          converged in "<<it_pf<<" iterations\n";
                    converged = true;
                    break;
                }
                else if (it_pf==maxIt-1)
                {
                    gsInfo<<"\t\tPF          did not converge!\n";
                    converged = false;
                    break;
                }

            } // NR phase field equation
            
            // Compute residual with updated damage (Dnew)

            auto D_pos = kappa * Heaviside_pos * reshape(PV,3,3) + 2*mu*reshape(PD,3,3);
            auto D_neg = kappa * Heaviside_neg * reshape(PV,3,3); // 3x3 matrix

            auto D_split_dmg = g_d * D_pos + D_neg;// material matrix with energy split (linear elasticity)
            K_u = igrad(w,G) * D_split_dmg * igrad(w,G).tr(); // stiffness matrix elastic problem
            
            // Compute new residual norm
            real_t res_stag =  (K_u * Unew - Q_u).norm();

            // Convergence check
            // if (stag == 0) res_d_norm0 = res_d;
            // else         res_d_norm = res_d;
            gsInfo<<"\t\tSTAG  iter   "<<stag<<": res_stag = "<<"\n";

            if (stag>0 && res_stag < tol_stag)
            {
                gsInfo<<"\t\tSTAG         converged in "<<stag<<" iterations\n";
                converged = true;
                break;
            }
            else if (stag==maxIt-1)
            {
                gsInfo<<"\t\tSTAG          did not converge!\n";
                converged = false;
                break;
            }
            // Check staggered convergence

        } // staggered iteration

        // Dnew = Dold;
        // Unew = Uold + delta_U; // update displacements for next step!

    } // step iteration

    if (plot)
    {
        collection.save();
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";

    return EXIT_SUCCESS;

} // end main
    
