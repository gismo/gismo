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
    ./bin/cahn-hilliard_example --plot -N 80 --clamped
    
    Run a simple Cahn-Hilliard example with an analytical initial condition "0.1 * cos(2*pi*x) * cos(2*pi*y)" (Nitsche) 
    ./bin/cahn-hilliard_example --plot -N 80 

    Run a simple Cahn-Hilliard example with an analytical initial condition "0.1 * cos(2*pi*x) * cos(2*pi*y)" until equilibrium (Nitsche)
    ./bin/cahn-hilliard_example --plot -N 1000  

    Mobility example (Gomez et al., 2014)
    ./bin/cahn-hilliard_example_mobility_works --plot -N 57   

    Mobility example (Random initial condition -- spinoidal decomposition)
    ./bin/cahn-hilliard_example_mobility_works --plot -N 100 --initial   
    
    
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
    real_t M0 = 1;
    // real_t dt = 1e-7;
    real_t dt = 1e-3;
    index_t maxSteps = 10;

    real_t eps_penalty = 1e4 *lambda;   // not sure about the value!


    bool output = true;


    index_t plotmod = 1;

    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 1;
    index_t numElevate = 1;
    bool last = false;
    bool nitsche = false;
    real_t mean = 0.0;
    // std::string fn("pde/ch_bvp_square.xml");

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
    index_t n = 40;
    index_t degree = 2;
    real_t hmax = 1.0/(n-degree);

    // 1. construction of a knot vector for each direction
    // n - degree - 1 interior knots
    // degree + 1 multiplicity at the ends
    // In total: n-p elements!
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
    bc.setGeoMap(mp);    

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
    gsMatrix<> dCnew,dCalpha,dCold, dCupdate, temp_conc;

    solution c = A.getSolution(w, Calpha); // C
    solution dc = A.getSolution(w, dCalpha); // \dot{C}

    // Mobility
    // auto M_c  = 1.0 + 0.0*c.val();
    // auto dM_c = 0.0;
    
    // Mobility and derivative (with respect to c)
    // For H. Gomez initial condition: ./bin/cahn-hilliard_mobility_example -N 123
    // auto M_c  = abs(M0 * c.val()*(1.0 - c.val()));
    // auto dM_c = M0 * (1.0 - 2.0*c.val()); // first derivative of M with respect to c
    // auto ddM_c = M0 * (- 2.0); // second derivative of M with respect to c
    
    // auto M_c  = ternary(1-abs(c.val()),1 + 0*c.val(), 1 + 0*c.val()); //abs(1.0 - (c*c).val());
    // // auto dM_c = (-2.0*c.val()); // first derivative of M with respect to c
    // // auto ddM_c = -2.0; // second derivative of M with respect to c
    // auto dM_c = ternary(1-abs(c.val()), 0*c.val(), 0*c.val()); // first derivative of M with respect to c
    // auto ddM_c = ternary(1-abs(c.val()), 0*c.val(), 0*c.val()); // second derivative of M with respect to c
    // real_t gamma_const = 0.87;// between 0 and 1
    
    // real_t gamma_const = 0.87;// between 0 and 1
    // auto M_c  = ternary(1-abs(c.val()), M0*(1 - pow(gamma_const,2)*(c*c).val()), 0*c.val()); //abs(1.0 - (c*c).val());
    // // auto dM_c = (-2.0*c.val()); // first derivative of M with respect to c
    // // auto ddM_c = -2.0; // second derivative of M with respect to c
    // auto dM_c = ternary(1-abs(c.val()), M0*(1 - pow(gamma_const,2)*2*c.val()), 0*c.val()); // first derivative of M with respect to c
    // auto ddM_c = ternary(1-abs(c.val()), M0*(-2*pow(gamma_const,2)) + 0*c.val(), 0*c.val()); // second derivative of M with respect to c

    // auto M_c  = ternary(1-abs(c.val()), M0*(1 - pow(gamma_const,2)*(c*c).val()), 0*c.val()); //abs(1.0 - (c*c).val());
    // // auto dM_c = (-2.0*c.val()); // first derivative of M with respect to c
    // // auto ddM_c = -2.0; // second derivative of M with respect to c
    // auto dM_c = ternary(1-abs(c.val()), M0*(1 - pow(gamma_const,2)*2*c.val()), 0*c.val()); // first derivative of M with respect to c
    // auto ddM_c = ternary(1-abs(c.val()), M0*(-2*pow(gamma_const,2)) + 0*c.val(), 0*c.val()); // second derivative of M with respect to c

    // auto M_c  = ternary(1-abs(c.val()), 1.0 + 0*c.val(), 0*c.val()); //abs(1.0 - (c*c).val());
    // // auto dM_c = (-2.0*c.val()); // first derivative of M with respect to c
    // // auto ddM_c = -2.0; // second derivative of M with respect to c
    // auto dM_c = ternary(1-abs(c.val()),0*c.val(),  0*c.val()); // first derivative of M with respect to c
    // auto ddM_c = ternary(1-abs(c.val()),0*c.val(),  0*c.val()); // second derivative of M with respect to c


    //auto M_c  = M0 * c.val()*(1.0 - c.val());
    // auto dM_c = M0 * (1.0 - 2.0*c.val()); // first derivative of M with respect to c
    // auto ddM_c = M0 * (- 2.0); // second derivative of M with respect to c

    // auto M_c  = M0 *(1.0 - (c*c).val());
    // auto dM_c = M0 * (1.0 - 2.0*c.val()); // first derivative of M with respect to c
    // auto ddM_c = M0 * (- 2.0); // second derivative of M with respect to c


    // auto M_c  = 1 + 0*c.val();
    // auto dM_c = 0 + 0*c.val(); // first derivative of M with respect to c
    // auto ddM_c = 0 + 0*c.val(); // second derivative of M with respect to c
    
    
    // auto M_c  = ; //abs(1.0 - (c*c).val());
    // // auto dM_c = (-2.0*c.val()); // first derivative of M with respect to c
    // // auto ddM_c = -2.0; // second derivative of M with respect to c
    // auto dM_c = ternary(1-abs(c.val()),0*c.val(),  0*c.val()); // first derivative of M with respect to c
    // auto ddM_c = ternary(1-abs(c.val()),0*c.val(),  0*c.val()); // second derivative of M with respect to c


    // Modified mobility
    // For H. Gomez initial condition:./bin/cahn-hilliard_mobility_example -N 60       
    // auto M_c  = abs(1.0 - (c*c).val());
    // auto M_c = if(abs(c.val())<1,1-pow(c,2),0);
    // auto M_c  = abs(1.0 - (c*c).val());

    auto M_c  = ternary(1-abs(c.val()), 1.0-(c*c).val(), 0*c.val()); //abs(1.0 - (c*c).val());
    // auto dM_c = (-2.0*c.val()); // first derivative of M with respect to c
    // auto ddM_c = -2.0; // second derivative of M with respect to c
    auto dM_c = ternary(1-abs(c.val()), -2.0*c.val(), 0*c.val()); // first derivative of M with respect to c
    auto ddM_c = ternary(1-abs(c.val()), -2.0 + 0*c.val(), 0*c.val()); // second derivative of M with respect to c

    // auto M_c  = (1.0 - (c*c).val());
    // auto dM_c = (-2.0*c.val()); // first derivative of M with respect to c
    // auto ddM_c = -2.0; // second derivative of M with respect to c

    // auto M_c  = M0 * c.val();
    // auto dM_c = M0 + 0. * c.val(); // first derivative of M with respect to c
    // auto ddM_c = 0.+ 0.* c.val(); // second derivative of M with respect to c

    // auto M_c  = M0 * c.val()*(1.0 - c.val());
    // auto dM_c = M0 * igrad(c,G) - 2.0 * M0 * igrad(c,G); 

    
    if (nitsche) 
    {   
        for ( gsMultiPatch<>::const_biterator
            bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
            {
                bc.addCondition( *bit, condition_type::neumann, nullptr); 
            }
        // A.initMatrix();
        // A.assembleBdr(bc.get("Neumann"), - lambda * M_c.val() * igrad(w,G) *  nv(G)  * ilapl(w,G).tr() + // consistency term
        //                                 (eps_penalty * hmax) * (igrad(w,G) * nv(G).normalized()) * (igrad(w,G) * nv(G)).tr()  // penalty (stabilizing) term
        //                                 - lambda * M_c.val() * ilapl(w,G) * (igrad(w,G)  * nv(G)).tr()); // symmetry term
        // K_nitsche = A.giveMatrix(); // .giveMatrix() moves the matrix A into K_nitche (avoids having two matrices A and K_nitsche)
        // gsInfo<<K_nitsche.toDense()<<"\n";
    }
    else
    {
        for ( gsMultiPatch<>::const_biterator
            bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
            {
                bc.addCondition( *bit, condition_type::clamped,0);
            }

    }

    // auto mu_c = 1.0 / (2.0*theta) * (c / (1.0-c).val()).log() + 1 - 2*c;
    // auto dmu_c= 1.0 / (2.0*theta) * igrad(c,G) / (c - c*c).val() - 2.0 * igrad(c,G);
    
    auto mu_c = pow(c,3).val() - c.val();
    auto dmu_c= igrad(c,G) * (- 1.0 + 3.0 * (c*c).val());

    // auto mu_c= -c.val() * (1.0 - (c*c).val());
    // auto dmu_c = -igrad(c,G) * (1.0 - (c*c).val()) - c.val() * (1.0 - 2.0 * c.val() * igrad(c,G));
    // auto M_c  = M0 * c * (1.0-c.val());
    // auto dM_c = M0 * igrad(c,G) - 2.0 * M0 * igrad(c,G);

    // Derivatives of the double well potential
    auto f_2 = - 1.0 + 3.0 * (c*c).val();
    auto f_3 = 6*c.val();

    // auto residual = w*dc + M_c.val()*igrad(w,G)*dmu_c.tr() +
    //                 lambda*ilapl(c,G).val()*igrad(w,G)*dM_c.tr() + M_c.val() * ilapl(w,G)*lambda*ilapl(c,G);
    // auto residual = w*dc + M_c.val()*igrad(w,G)*dmu_c.tr() +
    //                     lambda*ilapl(c,G).val()*igrad(w,G)*dM_c.tr() + M_c.val() * ilapl(w,G)*lambda*ilapl(c,G);
    // // auto residual = w*dc + // M
    //                 igrad(w,G)  * (- 1.0 + 3.0 * (c*c).val()) * igrad(c,G).tr() + // F_bar 
    //                 //igrad(w,G)  * (-1) * igrad(c,G).tr() + // F_bar 
    //                 // lambda*ilapl(c,G).val()*igrad(w,G)*dM_c.tr() + // term gradient mobility!
    //                 ilapl(w,G)*lambda*ilapl(c,G).val(); // K_laplacian
    
    auto residual = w*dc + // M
                    (f_2 * M_c.val()) * igrad(w,G) * igrad(c,G).tr() + // F_bar
                    (lambda * M_c.val()) * ilapl(w,G) * ilapl(c,G).val() + // K_laplacian 
                    (lambda * dM_c.val())* igrad(w,G) * ilapl(c,G).val() * igrad(c,G).tr(); // mobility term!

    //! [Problem setup]

    // Define linear solver (install SuperLUMT-devel)
#ifdef GISMO_WITH_SUPERLU
    gsSparseSolver<>::SuperLU solver;
#   else
    gsSparseSolver<>::LU solver;
#endif
    
    // Generalized-alpha method parameters
    real_t rho_inf = 0.5;
    real_t alpha_m = 0.5*(3-rho_inf) / (1+rho_inf);
    real_t alpha_f = 1 / (1+rho_inf);
    real_t gamma   = 0.5 + alpha_m - alpha_f;
    // time stepping options
    index_t maxIt = 100; // max NR iterations

    gsMatrix<> Q;
    gsSparseMatrix<> K, K_nitsche;
    

    // Legend:
    // C_old   = C_n
    // C       = C_n+1,i-1
    // C_alpha = C_{n+alpha_f,i-1}
    // dC

    gsInfo<<"Starting.."<<"\n";
    
    // Setup the space (compute Dirichlet BCs)
    w.setup(bc, dirichlet::l2Projection, 0);

    gsInfo<<"Initial condition.."<<"\n";
    
    if (random) 
    {
        // %%%%%%%%%%%%%%%%%%%%%%%% Random initial condition %%%%%%%%%%%%%%%%%%%%%%%%
        gsMatrix<> tmp = gsMatrix<>::Random(A.numDofs(),1);
        Cold = tmp.array()*0.1/2; //random uniform variable in [-0.05,0.05]
        Cold.array() += mean; // 0.45
    }
    else 
    {
        // %%%%%%%%%%%%%%%%%%%%%%%% Analytical intial condition %%%%%%%%%%%%%%%%%%%%%%%%
        gsFunctionExpr<> source = do3D ?
            gsFunctionExpr<>("0.1 * cos(2*pi*x) * cos(2*pi*y) * cos(2*pi*z)",3) :
            gsFunctionExpr<>("0.1 * cos(2*pi*x) * cos(2*pi*y)",2);
        gsMatrix<> tmp;
        Cold.setZero(A.numDofs(),1);
        real_t error = gsL2Projection<real_t>::projectFunction(dbasis,source,mp,tmp);  // 3rd arg has to be multipatch
        gsInfo << "L2 projection error "<<error<<"\n";
        for (index_t i = 0; i < dbasis.basis(0).size(); i++)
            if (w.mapper().is_free(i))
                Cold(w.mapper().index(i),0) = tmp(i,0);
    }
    
    Calpha = Cold;
    dCold.setZero(A.numDofs(),1);

    real_t Q0norm = 1, Qnorm = 10;
    real_t tol = 1e-4;

    gsParaviewCollection collection("ParaviewOutput/solution", &ev);
    collection.options().setSwitch("plotElements", true);
    collection.options().setInt("plotElements.resolution", 4);
    collection.options().setInt("numPoints",(do3D) ? 10000 : 1000);

    real_t dt_old = dt;
    real_t t_rho = 0.9;
    real_t t_err = 1;
    index_t lmax = 1;
    real_t TOL = 1e-3;
    std::vector<gsMatrix<>> Csols(2);

    real_t tmp_alpha_m = 1;
    real_t tmp_alpha_f = 1;
    real_t tmp_gamma   = 1;

    real_t time = 0;
    bool converged = false;

    A.initSystem(); // Initialize the system (outside the loops)

    for (index_t step = 0; step!=maxSteps; step++)
    {
        for (index_t dt_it = 0; dt_it != lmax; dt_it++)
        {
            gsInfo<<"Time step "<<step<<"/"<<maxSteps<<", iteration "<<dt_it<<": dt = "<<dt<<", [t_start,t_end] = ["<<time<<" , "<<time+dt<<"]"<<"\n";
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
                    
                    A.assemble(residual * meas(G));
                    Q = A.rhs();

                    if (nitsche) 
                    {      
                        A.assembleBdr(bc.get("Neumann"), - lambda * M_c.val() * igrad(w,G) *  nv(G)  * ilapl(w,G).tr() + // consistency term
                                                        (eps_penalty * hmax) * (igrad(w,G) * nv(G).normalized()) * (igrad(w,G) * nv(G)).tr() + // penalty (stabilizing) term
                                                        - lambda * M_c.val() * ilapl(w,G) * (igrad(w,G)  * nv(G)).tr()); // symmetry term
                        K_nitsche = A.giveMatrix(); // .giveMatrix() moves the matrix A into K_nitche (avoids having two matrices A and K_nitsche)

                        Q.noalias() += K_nitsche * Calpha; // add the residual term from Nitche (using the matrix )
                        // Old code lines:
                        // A.assembleBdr(bc.get("Neumann"), - igrad(w,G) * nv(G) * lambda * ilapl(c,G).val()); // consistency term
                        // A.assembleBdr(bc.get("Neumann"),  (igrad(w,G) * nv(G).normalized()) * hmax * eps_penalty * (igrad(c,G) * nv(G)) ); // penalty term
                        // A.assembleBdr(bc.get("Neumann"), - lambda * ilapl(w,G) * igrad(c,G) * nv(G)); // symmetry term
                    }
                    //if (step > 59) gsInfo<<M_c.val()<<"\n";
                    
                    if (it == 0) Q0norm = Q.norm();
                    else         Qnorm = Q.norm();

                    gsInfo<<"\t\tNR iter   "<<it<<": res = "<<Qnorm/Q0norm<<"\n";
                    
                    if (it>0 && Qnorm/Q0norm < tol)
                    {
                        gsInfo<<"\t\t"<<method<<"converged in "<<it<<" iterations\n";
                        converged = true;
                        break;
                    }
                    else if (it==maxIt-1)
                    {
                        gsInfo<<"\t\t"<<method<<"did not converge!\n";
                        converged = false;
                        break;
                    }

                    // // gsInfo<<"Assembly K_m\n";
                    // A.assembleJacobian( residual * meas(G), dc );
                    // K_m = tmp_alpha_m * A.matrix(); 

                    // // gsInfo<<"Assembly K_f\n";
                    // A.assembleJacobian( residual * meas(G), c );
                    // K_f = tmp_alpha_f * tmp_gamma * dt * A.matrix();

                    // A.assemble(M_c.val() * igrad(w,G).tr() * igrad(w,G) * (- 1.0 + 3.0 * (c*c).val()) +
                    //             lambda*ilapl(w,G).tr()*igrad(w,G)*dM_c.tr() +
                    //             M_c.val() * ilapl(w,G)*lambda*ilapl(w,G).tr());

                    // A.assemble(meas(G) * (igrad(w,G) * (- 1.0 + 3.0 * (c*c).val())* igrad(w,G).tr()  + // K_f1
                    //                     igrad(w,G) * ((6.0 * c.val()) * igrad(c,G).tr() * w.tr()) + // K_f2
                    //                     // lambda * igrad(w,G)*dM_c.tr()*ilapl(w,G).tr()   +  // K_mobility
                    //                     lambda * ilapl(w,G) * ilapl(w,G).tr())); // K_laplacian

                    //%% Assembly of the tangent stiffness matrix (K_m and K_f simultaneously) %%
                    A.clearMatrix(); // Resets to zero the values of the already allocated to matrix (LHS)
                    // A.assemble(meas(G) * (w*w.tr()*tmp_alpha_m +// K_m
                    //                     (tmp_alpha_f * tmp_gamma * dt)* (M_c.val() * f_2 *igrad(w,G) * igrad(w,G).tr() + // K_f1
                    //                     (M_c.val() * f_3) * igrad(w,G) * igrad(c,G).tr() * w.tr() + // K_f2
                    //                     (M_c.val() * lambda) * ilapl(w,G) * ilapl(w,G).tr()))); // K_laplacian          
                    //                 // K_mobility    

                    // A.initMatrix();
                    A.assemble(meas(G) * (w*w.tr()*tmp_alpha_m +// K_m
                                        (tmp_alpha_f * tmp_gamma * dt)* ( M_c.val() * f_2 *igrad(w,G) * igrad(w,G).tr() + // K_f1
                                        (M_c.val() * f_3) * igrad(w,G) * igrad(c,G).tr() * w.tr() + // K_f2
                                        (M_c.val() * lambda) * ilapl(w,G) * ilapl(w,G).tr() + // K_laplacian          
                                        (lambda * dM_c.val() * ilapl(c,G).val()) * igrad(w,G) * igrad(w,G).tr() +  // K_mob_3      
                                        (lambda * ddM_c * ilapl(c,G).val()) * igrad(w,G) * igrad(c,G).tr() * w.tr() +  // K_mob_2    
                                        (lambda * dM_c.val()) * igrad(w,G) * igrad(c,G).tr() * ilapl(w,G).tr() )));  // K_mob_1   
                                        // lambda * igrad(w,G)*dM_c.tr()*ilapl(w,G).tr()   +  // K_mobility
                    
                    K = A.matrix(); 
                    
                    if (nitsche) 
                        K += (tmp_alpha_f * tmp_gamma * dt) * K_nitsche; // add the Nitsche term to the stiffness matrix

                    
                    // A.clearMatrix();
                    // A.assembleJacobian( residual * meas(G), dc );
                    // K_m = tmp_alpha_m * A.matrix();

                    // A.assembleJacobian( residual * meas(G), c );
                    // K_f = tmp_alpha_f * tmp_delta * dt * A.matrix();    

                    // K = K_m + K_f + (tmp_alpha_f * tmp_gamma * dt) * K_nitsche;




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

                    // Check that all values from Cnew are within the concentration bounds
                    // gsMatrix<> tmp = gsMatrix<>::setZero(A.numDofs(),1);
                    // temp_conc.setZero(A.numDofs(),1);

                    // temp_conc = Cnew; // Cnew to a vector

                    // Cnew_vec within the range [-1, 1]
                    // for (int i = 0; i < Cnew.size(); ++i) 
                    // {
                    //     if (Cnew(i) > 1.0) 
                    //     {
                    //         Cnew(i) = 1.0;
                    //     } 
                    //     else if (Cnew(i) < -1.0) 
                    //     {
                    //         Cnew(i) = -1.0;
                    //     }
                    // }

                    // Convert the vector back to gsMatrix
                    // Cnew = temp_conc;
                }
                if (!converged)
                    break;

                // %% Switch to generalized-alpha parameters (k=1)
                tmp_alpha_m = alpha_m;
                tmp_alpha_f = alpha_f;
                tmp_gamma = gamma;

                // %% For time step adaptivity %%
                // Csols[k] = Cnew; // k=0: BE, k=1: alpha
            }

            // %%%%%%%%%% For time step adaptivity %%%%%%%%%%
            // if (converged)
            // {
            //     t_err = (Csols[0] - Csols[1]).norm() / (Csols[1]).norm();
            //     dt_old = dt;
            //     dt *= t_rho * math::sqrt(TOL / t_err);
            //     if (t_err < TOL)
            //         break;
            // }
            // else
            // {
            //     dt_old = dt;
            //     dt *= t_rho;
            // }
        }

        time += dt_old;
        Cold = Cnew;
        dCold = dCnew;

        //! [Export visualization in ParaView]
        if (plot && step % plotmod==0)
        {
            Calpha = Cnew;
            // collection.newTimeStep(&mp);
            collection.newTimeStep(&mp);
            collection.addField(c,"numerical solution");
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
