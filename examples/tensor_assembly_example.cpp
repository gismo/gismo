/** @file assembly_example.cpp

    @brief Tutorial on how to use the gsAssembler class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

# include <gismo.h>

# include <gsTensor/gsTensorFunction.h>
# include <gsTensor/gsTensorAssembler.h>


using namespace gismo;

int main(int argc, char *argv[])
{
    short_t D = 2;
    index_t numRefine = 0,
        numElevate = 0,
        degree = 2,
        knots = 0;
    
    gsCmdLine cmd("Tutorial on assemblying a Poisson problem.");
    cmd.addInt("d", "dimension", "Dimension 2d/3d", D);
    cmd.addInt("p", "degree", "Degree of the initial tensor product basis", degree);
    cmd.addInt("k", "knots", "Number of knots of the initial tensor product basis", knots);
    cmd.addInt("u", "uniform", "Number of uniform refinement", numRefine);
    cmd.addInt("e", "elevate", "Number of degree elevation steps", numElevate);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]
    
    gsTensorFunction<> RR(2,1); // creates a constant function = 1
    // gsDebugVar(RR);

    std::vector<gsKnotVector<> > knot_vector(D, gsKnotVector<>(0, 1, knots, degree + 1));
    gsBasis<>::uPtr basis = gsBSplineBasis<>::create(knot_vector);

    if (numElevate > 0)
        basis->degreeElevate(numElevate);
    for (int i = 0; i < numRefine; ++i)
        basis->uniformRefine();

    gsTensorAssembler<> ta;
    ta.compute(*basis, STIFFNESS, RR); // DIFFUSION ?
    // ta.compute(*basis, MASS, RR); // DIFFUSION ? bug for this line?

    gsFiberMatrix<real_t> mat = ta.kronecker().toFiberMatrix();

    // //-------- START TEST MATRIX-FREE SOLVER
    // auto & kp = ta.kronecker();
    // auto kpop = memory::make_shared_not_owned( (gsLinearOperator<>*)(&kp));
    // gsMatrix<> vec(kp.cols(),1);
    // vec.setRandom();
    // gsInfo << "Random vector :"<< vec.transpose() <<"\n";
    // gsMatrix<> prod;
    // kp.apply(vec, prod); // computes product kp*vec;
    // gsInfo << ( prod.transpose() ) <<"\n";
    // gsInfo << ( (mat.toSparseMatrix() * vec).transpose() ) <<"\n";
    //
    // gsVector<> diag;
    // kp.diagonal_into(diag);
    // gsDiagonalOp<real_t> dop(diag);
    //
    // //gsConjugateGradient<> PCG( kpop /*, dop*/ ); // setup CG
    // gsBiCgStab<> PCG( kpop /*, dop*/ ); // this works for non-symmetric
    // gsMatrix<> x(kp.rows(),1);
    // x.setZero();
    // PCG.solve(prod, x);
    // gsInfo << "Solution :"<< x.transpose() <<"\n";
    // //------------ END TEST MATRIX-FREE SOLVER
    
    // gsInfo << "Matrix size: \n" << mat.rows() <<"\n";
    // gsInfo << "Matrix : \n" << mat.toSparseMatrix().toDense() <<"\n";

    // comparison with standard assembly
    gsExprAssembler<> A(1, 1);
    std::string fn = "planar/unitsquare.xml";
    gsFileData<> fd(fn);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    gsMultiPatch<real_t>::uPtr mp = gsReadFile<>(fn);
    gsInfo << " Got" << *mp << " \n";

    gsStopwatch timer_standard;
        timer_standard.restart();
    gsExprAssembler<>::geometryMap G = A.getMap(*mp);
    A.setIntegrationDomain(basis->domain());
    auto u = A.getSpace(*basis);
    //auto v = A.getTestSpace(*basis, 1);
    //auto w = A.getCoeff(RR);
    A.initMatrix();
    A.assemble( igrad(u, G) * igrad(u, G).tr() * meas(G) );
    // gsInfo << "Standard assembly time: " << timer_standard.stop() << " seconds\n";
    real_t assembly_time_standard = timer_standard.stop();
    gsInfo << "Standard assembly time: " << assembly_time_standard << " seconds\n";

    // gsMatrix<> mat_std = A.giveMatrix().toDense();
    // gsInfo << "Standard assembly matrix: \n" << mat_std <<"\n";
    // gsInfo << "Difference: \n" << (mat.toSparseMatrix().toDense() - mat_std).cwiseAbs().maxCoeff() <<"\n";

    // // rank-1 density with sinusoidal variation
    // std::vector<gsMatrix<real_t>> skeleton(2);
    //
    index_t n1 = basis->component(0).size();
    index_t n2 = basis->component(1).size();
    //
    // skeleton[0].resize(n1,1);
    // skeleton[1].resize(n2,1);
    //
    // --- compute Greville points ---
    auto grevillePoints = [](const gsBSplineBasis<real_t>& b)
    {
        const gsKnotVector<real_t>& kv = b.knots();
        const int p = b.degree();
        const index_t n = b.size();

        gsVector<real_t> g(n);

        for (index_t i = 0; i < n; ++i)
        {
            real_t s = 0;
            for (int m = 1; m <= p; ++m)
                s += kv[i + m];
            g[i] = s / p;
        }
        return g;
    };

    const gsBSplineBasis<real_t>& b1 =
        dynamic_cast<const gsBSplineBasis<real_t>&>(basis->component(0));
    const gsBSplineBasis<real_t>& b2 =
        dynamic_cast<const gsBSplineBasis<real_t>&>(basis->component(1));

    gsVector<real_t> gx = grevillePoints(b1);
    gsVector<real_t> gy = grevillePoints(b2);
    //
    // // --- build rho ---
    // for (index_t i = 0; i < n1; ++i)
    // {
    //     skeleton[0](i,0) = 1.0 + 0.8 * std::sin(6.0 * M_PI * gx[i]);
    // }
    // for (index_t j = 0; j < n2; ++j)
    // {
    //    skeleton[1](j,0) = 1.0 + 0.8 * std::cos(6.0 * M_PI * gy[j]);
    // }
    // // skeleton[1].col(0).setOnes();
    //
    // gsTensorFunction<> rho(*basis, skeleton);

    ////////////// test tensor approximation //////////////
    // gsTensorBSplineBasis<2,real_t> prBasis(dynamic_cast<gsTensorBSplineBasis<2,real_t>>(*basis));
    const gsTensorBSplineBasis<2,real_t>& prBasis = dynamic_cast<const gsTensorBSplineBasis<2,real_t>&>(*basis);

    gsInfo << "Interpolation basis: " << prBasis << "\n";

    const gsMatrix<real_t> pts = prBasis.anchors();

    // rho(x,y) = 1 + 0.3 sin(2 pi x) cos(2 pi y)
    // gsFunctionExpr<real_t> rho_density(
    //     "1.0 + 0.8*sin(6*pi*x)*cos(8*pi*y)",
    //     2
    // );
    // gsFunctionExpr<real_t> rho_density(
    //     "1.0 + 0.9*(1.0 + tanh(50*(0.4 - sqrt((x-0.5)^2 + (y-0.5)^2))))",
    //     2
    // );
    // example circle interface
    gsFunctionExpr<real_t> rho_density(
        "1.0 + 50.0*exp(-100*pow(sqrt((x-0.5)^2 + (y-0.5)^2) - 0.3, 2))",
        2
    );

    gsMatrix<real_t> vals = rho_density.eval(pts);

    auto Afap = prBasis.interpolateAtAnchors(vals);

    real_t eps = 1e-3;
    index_t maxrank = 10;

    gsTensorFunction<real_t> rho(*Afap, eps, maxrank);
    ////////////////////////////////////







    gsStopwatch timer_diffusion;
    timer_diffusion.restart();
    gsTensorAssembler<> ta_diffusion;
    ta_diffusion.compute(*basis, DIFFUSION, rho);
    real_t assembly_time_tensor = timer_diffusion.stop();
    gsInfo << "Tensor assembly time: " << assembly_time_tensor << " seconds\n";
    gsInfo << "Speedup: " << assembly_time_standard / assembly_time_tensor << "x\n";

    gsFiberMatrix<real_t> mat_diffusion = ta_diffusion.kronecker().toFiberMatrix();

    // gsInfo << "Diffusion matrix size:\n" << mat_diffusion.rows() << "\n";
    // gsInfo << "Diffusion matrix:\n" << mat_diffusion.toSparseMatrix().toDense() << "\n";

    // add boundary condition and solve
    // ------------------------------------------------------------
    // Apply Dirichlet BC directly on mat_diffusion and solve
    // ------------------------------------------------------------

    gsSparseMatrix<real_t> K = mat_diffusion.toSparseMatrix();
    const index_t N = K.rows();

    gsMatrix<real_t> rhs_x = gsMatrix<real_t>::Zero(N,1);
    gsMatrix<real_t> rhs_y = gsMatrix<real_t>::Zero(N,1);

    // tensor-product basis sizes
    // index_t n1 = basis->component(0).size();
    // index_t n2 = basis->component(1).size();

    GISMO_ASSERT(n1 * n2 == N, "Unexpected tensor-product size mismatch.");

    // ------------------------------------------------------------
    // Build Greville abscissae in each direction
    // ------------------------------------------------------------
    // gsVector<real_t> gx(n1), gy(n2);
    //
    // for (index_t i = 0; i < n1; ++i)
    //     gx[i] = basis->component(0).anchor(i).value();
    // for (index_t j = 0; j < n2; ++j)
    //     gy[j] = basis->component(1).anchor(j).value();

    // ------------------------------------------------------------
    // Collect boundary dofs and prescribed values
    // ------------------------------------------------------------
    std::vector<index_t> bdofs;
    std::vector<real_t> vals_x, vals_y;
    std::vector<bool> isBoundary(N, false);

    auto add_bc = [&](index_t i, index_t j)
    {
        index_t idx = i + n1 * j; // assuming first direction varies fastest
        if (isBoundary[idx]) return;
        isBoundary[idx] = true;

        bdofs.push_back(idx);
        vals_x.push_back(gx[i]); // x = xi
        vals_y.push_back(gy[j]); // y = eta
    };

    // left / right boundaries
    for (index_t j = 0; j < n2; ++j)
    {
        add_bc(0,      j);
        add_bc(n1 - 1, j);
    }

    // bottom / top boundaries
    for (index_t i = 0; i < n1; ++i)
    {
        add_bc(i, 0);
        add_bc(i, n2 - 1);
    }
    //
    // auto applyDirichletDense = [&](gsMatrix<real_t>& A,
    //                            gsMatrix<real_t>& rhs,
    //                            const std::vector<index_t>& dofs,
    //                            const std::vector<real_t>& vals)
    // {
    //     GISMO_ASSERT(dofs.size() == vals.size(), "Dirichlet size mismatch.");
    //
    //     for (size_t k = 0; k < dofs.size(); ++k)
    //     {
    //         index_t c = dofs[k];
    //         real_t g  = vals[k];
    //
    //         // rhs <- rhs - A(:,c) * g
    //         rhs -= A.col(c) * g;
    //
    //         // zero row and column
    //         A.row(c).setZero();
    //         A.col(c).setZero();
    //
    //         A(c,c) = 1.0;
    //         rhs(c,0) = g;
    //     }
    // };
    //
    // gsMatrix<real_t> Kx_dense = gsMatrix<real_t>(K);
    // gsMatrix<real_t> Ky_dense = gsMatrix<real_t>(K);
    //
    // applyDirichletDense(Kx_dense, rhs_x, bdofs, vals_x);
    // applyDirichletDense(Ky_dense, rhs_y, bdofs, vals_y);
    //
    // gsMatrix<real_t> sol_x = Kx_dense.fullPivLu().solve(rhs_x);
    // gsMatrix<real_t> sol_y = Ky_dense.fullPivLu().solve(rhs_y);

        auto solveWithDirichletSparse =
        [&](const gsSparseMatrix<real_t>& Kfull,
            const gsMatrix<real_t>& rhs_full,
            const std::vector<index_t>& bdofs,
            const std::vector<real_t>& vals) -> gsMatrix<real_t>
    {
        GISMO_ASSERT(bdofs.size() == vals.size(), "Dirichlet size mismatch.");

        const index_t N = Kfull.rows();
        GISMO_ASSERT(Kfull.cols() == N, "Matrix must be square.");
        GISMO_ASSERT(rhs_full.rows() == N && rhs_full.cols() == 1,
                     "rhs must be N x 1.");

        // Mark boundary dofs
        std::vector<bool> isBoundary(N, false);
        for (size_t k = 0; k < bdofs.size(); ++k)
        {
            GISMO_ASSERT(bdofs[k] >= 0 && bdofs[k] < N, "Boundary dof out of range.");
            isBoundary[bdofs[k]] = true;
        }

        // Build free dof list and old->new index map
        std::vector<index_t> freeDofs;
        freeDofs.reserve(N - bdofs.size());

        std::vector<index_t> oldToFree(N, -1);
        for (index_t i = 0; i < N; ++i)
        {
            if (!isBoundary[i])
            {
                oldToFree[i] = static_cast<index_t>(freeDofs.size());
                freeDofs.push_back(i);
            }
        }

        const index_t Nf = static_cast<index_t>(freeDofs.size());

        // Full prescribed vector ub
        gsMatrix<real_t> u = gsMatrix<real_t>::Zero(N, 1);
        for (size_t k = 0; k < bdofs.size(); ++k)
            u(bdofs[k], 0) = vals[k];

        // Reduced rhs: rf = f_f - K_fb * u_b
        gsMatrix<real_t> rhs_reduced(Nf, 1);
        rhs_reduced.setZero();

        for (index_t ii = 0; ii < Nf; ++ii)
        {
            const index_t iOld = freeDofs[ii];
            rhs_reduced(ii, 0) = rhs_full(iOld, 0);
        }

        rhs_reduced.noalias() -= (Kfull * u)(freeDofs, gsEigen::all);

        // Build Kff
        std::vector<gsEigen::Triplet<real_t>> trips;
        trips.reserve(Kfull.nonZeros());

        for (index_t col = 0; col < Kfull.outerSize(); ++col)
        {
            for (gsSparseMatrix<real_t>::InnerIterator it(Kfull, col); it; ++it)
            {
                const index_t iOld = it.row();
                const index_t jOld = it.col();

                if (!isBoundary[iOld] && !isBoundary[jOld])
                {
                    const index_t iNew = oldToFree[iOld];
                    const index_t jNew = oldToFree[jOld];
                    trips.emplace_back(iNew, jNew, it.value());
                }
            }
        }

        gsSparseMatrix<real_t> Kff(Nf, Nf);
        Kff.setFromTriplets(trips.begin(), trips.end());
        Kff.makeCompressed();

        // Solve reduced system
        gsEigen::SimplicialLDLT<gsSparseMatrix<real_t>> solver;
        solver.compute(Kff);
        GISMO_ASSERT(solver.info() == gsEigen::Success,
                     "Factorization failed for reduced sparse system.");

        gsMatrix<real_t> uf = solver.solve(rhs_reduced);
        GISMO_ASSERT(solver.info() == gsEigen::Success,
                     "Solve failed for reduced sparse system.");

        // Scatter back to full vector
        for (index_t ii = 0; ii < Nf; ++ii)
            u(freeDofs[ii], 0) = uf(ii, 0);

        return u;
    };

    gsMatrix<real_t> rhs0 = gsMatrix<real_t>::Zero(K.rows(), 1);

    gsMatrix<real_t> sol_x = solveWithDirichletSparse(K, rhs0, bdofs, vals_x);
    gsMatrix<real_t> sol_y = solveWithDirichletSparse(K, rhs0, bdofs, vals_y);

    // for (size_t k = 0; k < bdofs.size(); ++k)
    // {
    //     gsInfo << "Boundary dof: " << bdofs[k]
    //            << ", prescribed x = " << vals_x[k]
    //            << ", prescribed y = " << vals_y[k] << "\n";
    // }

    // recover b-spline solution
    // ------------------------------------------------------------
    // Reconstruct B-spline geometry from solved x/y coefficients
    // ------------------------------------------------------------
    GISMO_ASSERT(sol_x.rows() == N && sol_x.cols() == 1, "sol_x has wrong size.");
    GISMO_ASSERT(sol_y.rows() == N && sol_y.cols() == 1, "sol_y has wrong size.");

    gsMatrix<real_t> coefs(N, 2);
    coefs.col(0) = sol_x;
    coefs.col(1) = sol_y;

    // Recover geometry: X(xi,eta) = (x(xi,eta), y(xi,eta))
    gsGeometry<real_t>::uPtr geom = basis->makeGeometry(coefs);
    // gsInfo << "Recovered geometry:\n" << *geom << "\n";

    // Put geometry into a multipatch for Paraview output
    gsMultiPatch<real_t> mp_out;
    mp_out.addPatch(*geom);

    // Write Paraview
    gsWriteParaview(mp_out, "weighted_laplace_solution", 1000);

    return 0;
}

