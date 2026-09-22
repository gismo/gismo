/** @file fischerkolmogrov_prostate_example.cpp

    @brief Tutorial on how to use G+Smo to solve a Fischer-Kolmogorov problem.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    mkdir -p 00_output_FK_cube
    ./build/bin/fischerkolmogrov_prostate_example \
        -f pde/fisher_kolmogorov_cube_bvp.xml -r 4 -N 10 -t 0.5 \
        -D 1e-3 -p 0.5 -R 0.15 -s 0.075 --Nth 0.15 \
        -o 00_output_FK_cube --plot

    Author(s): L. Venta-Viñuela
*/
#include <gismo.h>

using namespace gismo;


// =============================================================================
// FISHER-KOLMOGOROV MODEL PARAMETERS (user-settable, for literature calibration)
// =============================================================================
// The values below are only the fall-back defaults, used when the label is not
// found in the xml file (see main()). They are chosen for a domain of unit size;
// D and rho carry the units of the geometry, so for e.g. a prostate measured in
// mm and time in days they have to be re-calibrated ([mm^2/day] and [1/day]).
template <class T>
struct fkParams
{
    T D     = 1.0e-3; // tumor cell diffusivity          [length^2/time]
    T rho   = 0.5;    // net proliferation rate          [1/time]
    T R     = 0.15;   // radius of the initial tumor     [length]
    T sigma = 0.05;   // interface width of the IC       [length]
    T Nth   = 0.15;   // threshold defining Omega_T      [-]
    // Centre of the initial tumor (physical coordinates).
    // NaN means: use the centre of the bounding box of the domain.
    T xc = std::numeric_limits<T>::quiet_NaN();
    T yc = std::numeric_limits<T>::quiet_NaN();
    T zc = std::numeric_limits<T>::quiet_NaN();
};

// Smooth spherical tumor seed:  N(x,0) = 0.5*(1 - tanh((||x - x_tumor|| - R)/sigma))
template <short_t dim, class T>
gsFunctionExpr<T> tumorSeed(const fkParams<T> & par, const gsMultiPatch<T> & mp)
{
    // Default seed location: centre of the bounding box of the domain
    gsMatrix<T> bbox;
    mp.boundingBox(bbox);
    gsVector<T,3> c;
    c.setZero();
    for (short_t d = 0; d != dim; d++)
        c[d] = 0.5 * (bbox(d,0) + bbox(d,1));

    if (!math::isnan(par.xc))            c[0] = par.xc;
    if (!math::isnan(par.yc))            c[1] = par.yc;
    if (dim == 3 && !math::isnan(par.zc)) c[2] = par.zc;

    std::ostringstream expr;
    expr.precision(16);
    expr << "0.5*(1-tanh((sqrt((x-(" << c[0] << "))^2+(y-(" << c[1] << "))^2";
    if (dim == 3)
        expr << "+(z-(" << c[2] << "))^2";
    expr << ")-(" << par.R << "))/(" << par.sigma << ")))";

    gsInfo << "Initial condition: N(x,0) = " << expr.str() << "\n";
    return gsFunctionExpr<T>(expr.str(), dim);
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
           index_t & numRefine,
           const fkParams<T> & par,
           std::string out,
           bool plot)
{
    // -------------------------------------------------------------------------
    // 1. Basis & Problem Setup (Uniform Tensor Product Basis, No Adaptivity)
    // -------------------------------------------------------------------------
    gsMultiBasis<> dbasis(mp, true);
    for (index_t r = 0; r != numRefine; ++r)
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
    
    // Parameters (fixed model constants, set from the command line)
    real_t D   = par.D;   // tumor cell diffusivity
    real_t rho = par.rho; // net proliferation rate
    gsInfo << "Fisher-Kolmogorov parameters: D = " << D << ", rho = " << rho
           << ", R = " << par.R << ", sigma = " << par.sigma << ", N_th = " << par.Nth << "\n";
    gsInfo << "Discretization: numRefine = " << numRefine << ", DoFs = " << dbasis.size()
           << ", dt = " << dt << ", steps = " << maxSteps << "\n";

    // Initialization biomarkers
    real_t V_T, V_N, V_P, N_dash, N_T, A_p;
    real_t int_N, int_H; // \int_Omega N dOmega  and  \int_Omega H(N-N_th) dOmega

    std::ofstream csvFile;
    csvFile.open(out+"/output_data.csv");
    csvFile << "TimeStep,Time,NumDOFs,int_N,int_H,V_P,N_dash,N_T,A_p\n";
    
    real_t N_hat_th = par.Nth;  // Cell density threshold defining Omega_T (tumor)
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

    // Domain residual of  dN/dt = div(D grad N) + rho*N*(1-N):
    //   \int w*dN + \int D*grad(w).grad(N) - \int rho*w*N*(1-N) = 0
    auto residual = w * dN +
                    D * igrad(w, G) * gradN.tr() +
                    (-rho) * w * N.val() * (1.0 - N.val());

    w.setup(bc, dirichlet::l2Projection, 0);
    A.initSystem();

    gsSparseSolver<>::uPtr solver = gsSparseSolver<>::get("SimplicialLDLT");
    gsSparseMatrix<> K, K_const, K_nitsche;

    // -------------------------------------------------------------------------
    // 4. Initial Condition Setup (smooth hyperbolic-tangent tumor seed)
    // -------------------------------------------------------------------------
    gsFunctionExpr<T> N0 = tumorSeed<dim,T>(par, mp);

    // L2-projection of the (physical-space) seed onto the solution basis
    gsMatrix<T> N0coefs;
    gsInfo << "IC L2-projection error: "
           << gsL2Projection<T>::project(dbasis, mp, N0, N0coefs) << "\n";

    mp_nold.addPatch(dbasis.basis(0).makeGeometry(N0coefs));
    mp_dnold.addPatch(dbasis.basis(0).makeGeometry(gsMatrix<T>::Zero(dbasis.basis(0).size(),1)));

    Nold.setZero(A.numDofs(), 1);
    dNold.setZero(A.numDofs(), 1);
    for (index_t i = 0; i != dbasis.basis(0).size(); i++)
        if (w.mapper().is_free(i))
            Nold(w.mapper().index(i),0) = N0coefs(i,0);

    // -------------------------------------------------------------------------
    // 5. Pre-assemble Constant Terms Matrix
    // -------------------------------------------------------------------------
    A.assemble(meas(G) * (w * w.tr() * alpha_m +
               (alpha_f * gamma * dt) * (D * igrad(w, G) * igrad(w, G).tr())));
    K_const = A.giveMatrix();

    // -------------------------------------------------------------------------
    // 6. ParaView output of the cell density N (one time step per frame)
    // -------------------------------------------------------------------------
    gsExprEvaluator<> ev(A);
    gsParaviewCollection collection(out + "/N", &ev);
    collection.options().setSwitch("plotElements", false); // the mesh is static, do not rewrite it
    collection.options().setInt("numPoints", (dim == 3) ? 50000 : 5000);
    if (plot)
    {   // initial condition; frame i of the collection = row i-1 of the csv file
        collection.newTimeStep(&mp);
        collection.addField(nold, "N");
        collection.saveTimeStep();
    }

    // -------------------------------------------------------------------------
    // 7. Time Integration & Newton-Raphson Solving Loop
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
            A.assemble(meas(G) * (-alpha_f * gamma * dt) *
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
        Nold  = Nnew;
        dNold = dNnew;

        if (plot)
        {
            collection.newTimeStep(&mp);
            collection.addField(nnew_sol, "N");
            collection.saveTimeStep();
        }
        
        // Indicator expression for tumor domain Omega_T (where N >= N_th) | N_th is 0.15
        auto tumor_indicator = 1.0 / (1.0 + (-(k_steep) * (nnew_sol.val() - N_hat_th)).exp());

        // ==============================================================
        // Quantities of interest of the forward simulation
        int_N = computeInt(dbasis, mp, nnew_sol.val(), A.options());  // \int_Omega N dOmega
        int_H = computeInt(dbasis, mp, tumor_indicator, A.options()); // \int_Omega H(N-N_th) dOmega
        gsInfo << "  int_Omega N dOmega          = " << int_N << "\n"
               << "  int_Omega H(N-N_th) dOmega  = " << int_H << "\n";

        // Compute biomarkers
        V_T = int_H; // Tumor volume
        auto VN_field = nnew_sol.val() * tumor_indicator; // needed tumor_indicator because we integrate over Omega_T
        V_N = computeInt(dbasis, mp, VN_field, A.options()); // total tumor cell volume

        V_P = computeInt(dbasis, mp, 1.0, A.options()); // prostate volume

        N_dash = V_N / V_T; // mean normalized tumor cell density
        N_T = N_dash * V_T / V_P; // total tumor index
        auto A_p_field = rho * nnew_sol.val() * (1.0 - nnew_sol.val()) * tumor_indicator; //we integrate over Omega_T
        A_p = computeInt(dbasis,mp,A_p_field,A.options()); // mean proliferation activity of the tumor
        // ==============================================================

        // Write in a csv file
        csvFile << step << "," << (step+1)*dt << "," << dbasis.size() << ","
                << int_N << "," << int_H << ","
                << V_P << "," << N_dash << "," << N_T << "," << A_p << "\n";
        csvFile.flush();
    }
    csvFile.close();
    if (plot) collection.save();

}

// =============================================================================
// MAIN FUNCTION
// =============================================================================
int main(int argc, char *argv[])
{
    std::string fn("pde/fisher_kolmogorov_cube_bvp.xml");
    std::string out("output");

    // Pre-scan the command line for the setup file, so that the option lists it
    // contains can serve as the defaults of all the other arguments below.
    // Precedence: built-in default  <  xml file  <  command line.
    for (int i = 1; i+1 != argc; ++i)
    {
        const std::string arg(argv[i]);
        if (arg == "-f" || arg == "--file") fn = argv[i+1];
    }

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

    // Model parameters (xml id 3), the labels match the long command line names
    fkParams<real_t> par;
    par.D     = CHopt.askReal("diffusivity", par.D    );
    par.rho   = CHopt.askReal("rho",         par.rho  );
    par.R     = CHopt.askReal("radius",      par.R    );
    par.sigma = CHopt.askReal("sigma",       par.sigma);
    par.Nth   = CHopt.askReal("Nth",         par.Nth  );
    par.xc    = CHopt.askReal("xc",          par.xc   );
    par.yc    = CHopt.askReal("yc",          par.yc   );
    par.zc    = CHopt.askReal("zc",          par.zc   );

    // Time integration (xml id 4) and discretization (xml id 5)
    real_t  dt        = TIMEopt.askReal  ("dt",        0.5);
    index_t maxSteps  = TIMEopt.askInt   ("Nsteps",    10 );
    index_t numRefine = Aopt.askInt      ("numRefine", 4  );
    bool    plot      = Aopt.askSwitch   ("plot",   false );

    // Command-line interface
    gsCmdLine cmd("Fisher-Kolmogorov Generalized-Alpha Solver (forward problem)");
    cmd.addReal("t", "dt", "Time step size", dt);
    cmd.addInt("N", "Nsteps", "Number of time steps", maxSteps);
    cmd.addInt("r", "numRefine", "Number of uniform refinement steps", numRefine);
    cmd.addString("f", "file", "Input XML setup file", fn);
    cmd.addString( "o", "output", "Output directory", out);
    cmd.addSwitch("plot", "Write ParaView files of N for every time step", plot);
    cmd.addReal("D", "diffusivity", "Tumor cell diffusivity D [length^2/time]", par.D);
    cmd.addReal("p", "rho", "Net proliferation rate rho [1/time]", par.rho);
    cmd.addReal("R", "radius", "Radius R of the initial tumor [length]", par.R);
    cmd.addReal("s", "sigma", "Interface width sigma of the initial tumor [length]", par.sigma);
    cmd.addReal("", "Nth", "Cell density threshold N_th defining Omega_T", par.Nth);
    cmd.addReal("", "xc", "x-coordinate of the tumor centre (default: domain centre)", par.xc);
    cmd.addReal("", "yc", "y-coordinate of the tumor centre (default: domain centre)", par.yc);
    cmd.addReal("", "zc", "z-coordinate of the tumor centre (default: domain centre)", par.zc);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    // Dispatch solver based on domain dimension
    if (mp.geoDim() == 2)
        solve<2, real_t>(mp, bc, CHopt, TIMEopt, Aopt, dt, maxSteps, numRefine, par, out, plot);
    else if (mp.geoDim() == 3)
        solve<3, real_t>(mp, bc, CHopt, TIMEopt, Aopt, dt, maxSteps, numRefine, par, out, plot);
    else
        GISMO_ERROR("Only 2D and 3D geometries are supported.");

    return EXIT_SUCCESS;
}