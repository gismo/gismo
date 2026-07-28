#include <gismo.h>

using namespace gismo;

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

// =============================================================================
// SOLVER ROUTINE
// =============================================================================
template <short_t dim, class T>
void solve(gsMultiPatch<T> & mp,
           gsBoundaryConditions<T> & bc,
           gsOptionList & CHopt,
           gsOptionList & TIMEopt,
           gsOptionList & Aopt,
           real_t & dt,
           index_t & maxSteps,
           std::string ic_file,
           std::string out)
{
    // -------------------------------------------------------------------------
    // 1. Basis & Problem Setup (Uniform Tensor Product Basis, No Adaptivity)
    // -------------------------------------------------------------------------
    gsMultiBasis<> dbasis(mp, true);
    for (size_t p = 0; p != dbasis.nBases(); ++p)
    {
        dbasis.basis(p).uniformRefine();
    }

    real_t hmax    = dbasis.basis(0).getMaxCellLength();
    real_t lambda  = CHopt.askReal("lambda", 1.0 / (32.0 * pow(EIGEN_PI, 2)));
    real_t penalty = 1e4 * lambda;

    // -------------------------------------------------------------------------
    // 2. Time Integration Setup (Generalized-Alpha Method)
    // -------------------------------------------------------------------------
    real_t rho_inf = TIMEopt.askReal("rho_inf", 0.5);
    real_t alpha_m = 0.5 * (3.0 - rho_inf) / (1.0 + rho_inf);
    real_t alpha_f = 1.0 / (1.0 + rho_inf);
    real_t gamma   = 0.5 + alpha_m - alpha_f;
    
    // Parameters
    real_t D = 0.0; //tumor cell diffusivity  
    real_t rho = 0.0; //net proliferation rate
    
    // Initialization biomarkers
    real_t V_T, V_N, V_P, N_dash, N_T, A_p;

    std::ofstream csvFile;
    csvFile.open(out+"/output_data.csv");
    csvFile << "TimeStep,NumDOFs,V_P,N_dash,N_T,A_p\n";
    
    real_t N_hat_th = 0.15;  // Cell density threshold defining Omega_T
    real_t k_steep = 200.0;  // Steepness factor for Heaviside step function


    gsConstantFunction<> alpha_f_func(alpha_f, dim);
    gsConstantFunction<> alpha_m_func(alpha_m, dim);

    // -------------------------------------------------------------------------
    // 3. Assembler & Expression Setup
    // -------------------------------------------------------------------------
    gsExprAssembler<> A(1, 1);
    A.options().setSwitch("SameElement", Aopt.askSwitch("SameElement", true));
    A.setIntegrationElements(dbasis);

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    geometryMap G = A.getMap(mp);
    space       w = A.getSpace(dbasis);

    gsMatrix<> Nnew, dNnew, Nold, dNold;
    gsMatrix<> Nalpha, dNalpha, dNupdate, Q;

    gsMultiPatch<> mp_nold, mp_dnold;

    auto nold  = A.getCoeff(mp_nold);
    auto dnold = A.getCoeff(mp_dnold);

    solution nnew_sol  = A.getSolution(w, Nnew);
    solution dnnew_sol = A.getSolution(w, dNnew);

    auto af = A.getCoeff(alpha_f_func);
    auto am = A.getCoeff(alpha_m_func);

    // Interpolated state variables
    auto N    = nold  + af * (nnew_sol - nold); // N is the tumor cell density
    auto dN    = dnold + am * (dnnew_sol - dnold);
    auto gradN = igrad(nold, G) + af * (igrad(nnew_sol, G) - igrad(nold, G));

    // Domain residual
    auto residual = w * dN +
                    D * igrad(w, G) * gradN.tr() +
                    rho * w * N.val() * (1.0 - N.val());

    w.setup(bc, dirichlet::l2Projection, 0);
    A.initSystem();

    gsSparseSolver<>::uPtr solver = gsSparseSolver<>::get("SimplicialLDLT");
    gsSparseMatrix<> K, K_const, K_nitsche;

    // -------------------------------------------------------------------------
    // 4. Initial Condition Setup (Random IC from XML File)
    // -------------------------------------------------------------------------
    gsFileData<> fd;
    fd.read(ic_file);
    gsMultiPatch<> MP;
    fd.getId(0, MP);

    mp_nold.addPatch(MP.patch(0));
    Nold.setZero(dbasis.size(), 1);
    mp_dnold.addPatch(dbasis.basis(0).makeGeometry(Nold));

    // -------------------------------------------------------------------------
    // 5. Pre-assemble Constant Terms Matrix
    // -------------------------------------------------------------------------
    A.assemble(meas(G) * (w * w.tr() * alpha_m +
               (alpha_f * gamma * dt) * (D * igrad(w, G) * igrad(w, G).tr())));
    K_const = A.giveMatrix();

    // -------------------------------------------------------------------------
    // 6. Time Integration & Newton-Raphson Solving Loop
    // -------------------------------------------------------------------------
    index_t maxIt = 50;
    real_t tol    = TIMEopt.askReal("tol", 1e-4);

    for (index_t step = 0; step < maxSteps; ++step)
    {
        gsInfo << "Time Step " << step + 1 << " / " << maxSteps << "\n";

        Nnew  = Nold;
        dNnew = (gamma - 1.0) / gamma * dNold;

        for (index_t it = 0; it < maxIt; ++it)
        {
            Nalpha  = Nold  + alpha_f * (Nnew  - Nold);
            dNalpha = dNold + alpha_m * (dNnew - dNold);

            // Domain Residual Assembly
            A.initMatrix();
            A.clearRhs();
            A.assemble(residual * meas(G));
            Q = A.rhs();

            if (it > 0 && Q.norm() < tol)
            {
                gsInfo << "  Converged in " << it << " iterations.\n";
                break;
            }

            // Tangent Matrix Assembly
            A.initMatrix();
            A.assemble(meas(G) * (alpha_f * gamma * dt) *
                                 (rho * w * (1.0 - 2.0 * N.val()) * w.tr()));

            K = A.giveMatrix();
            K += K_const;

            // Solve and Update
            solver->compute(K);
            dNupdate = solver->solve(-Q);

            dNnew += dNupdate;
            Nnew.noalias() += (gamma * dt) * dNupdate;
        }

        // Advance step
        nnew_sol.extract(mp_nold);
        dnnew_sol.extract(mp_dnold);
        
        // Indicator expression for tumor domain Omega_T (where N >= N_th) | N_th is 0.15
        auto tumor_indicator = 1.0 / (1.0 + (-(k_steep) * (nnew_sol.val() - N_hat_th)).exp());

        // ==============================================================
        // Compute biomarkers
        V_T = computeInt(dbasis, mp, tumor_indicator, A.options());; // Tumor volume
        auto VN_field = nnew_sol.val() * tumor_indicator; // needed tumor_indicator because we integrate over Omega_T
        V_N = computeInt(dbasis, mp, VN_field, A.options()); // total tumor cell volume

        V_P = computeInt(dbasis, mp, 1.0, A.options()); // prostate volume

        N_dash = V_N / V_T; // mean normalized tumor cell density
        N_T = N_dash * V_T / V_P; // total tumor index
        auto A_p_field = rho * nnew_sol.val() * (1.0 - nnew_sol.val()) * tumor_indicator; //we integrate over Omega_T
        A_p = computeInt(dbasis,mp,A_p_field,A.options()); // mean proliferation activity of the tumor
        // ==============================================================

        // Write in a csv file
        csvFile << step << "," << dbasis.size() << "," << V_P << "," << N_dash << "," << N_T << "," << A_p << "\n";
        csvFile.flush();
    }
    csvFile.close();

}

// =============================================================================
// MAIN FUNCTION
// =============================================================================
int main(int argc, char *argv[])
{
    real_t dt        = 1e-3;
    index_t maxSteps = 10;
    std::string fn("pde/cahn_hilliard_bvp.xml");
    std::string ic_file("random_ic/ic_nucleation.xml");
    std::string out("output");

    // Command-line interface
    gsCmdLine cmd("Cahn-Hilliard Generalized-Alpha Solver (No Adaptivity)");
    cmd.addReal("t", "dt", "Time step size", dt);
    cmd.addInt("N", "Nsteps", "Number of time steps", maxSteps);
    cmd.addString("f", "file", "Input XML setup file", fn);
    cmd.addString("i", "icfile", "Initial Condition XML file", ic_file);
    cmd.addString( "o", "output", "Output directory", out);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    // Read input problem settings
    gsFileData<> fd(fn);
    gsInfo << "Loaded configuration file: " << fd.lastPath() << "\n";

    gsMultiPatch<> mp;
    fd.getId(0, mp); // Spatial geometry / domain

    gsBoundaryConditions<> bc;
    fd.getId(2, bc); // Boundary conditions
    bc.setGeoMap(mp);

    gsOptionList CHopt, TIMEopt, Aopt;
    fd.getId(3, CHopt);   // Physics options
    fd.getId(4, TIMEopt); // Integrator options
    fd.getId(5, Aopt);    // Assembler options

    // Dispatch solver based on domain dimension
    if (mp.geoDim() == 2)
        solve<2, real_t>(mp, bc, CHopt, TIMEopt, Aopt, dt, maxSteps, ic_file, out);
    else if (mp.geoDim() == 3)
        solve<3, real_t>(mp, bc, CHopt, TIMEopt, Aopt, dt, maxSteps, ic_file, out);
    else
        GISMO_ERROR("Only 2D and 3D geometries are supported.");

    return EXIT_SUCCESS;
}