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
    auto v = A.getTestSpace(*basis, 1);
    auto w = A.getCoeff(RR);
    A.initMatrix();
    A.assemble( igrad(u, G) * igrad(u, G).tr() * meas(G) );
    // gsInfo << "Standard assembly time: " << timer_standard.stop() << " seconds\n";
    real_t assembly_time_standard = timer_standard.stop();
    gsInfo << "Standard assembly time: " << assembly_time_standard << " seconds\n";

    gsMatrix<> mat_std = A.giveMatrix().toDense();
    // gsInfo << "Standard assembly matrix: \n" << mat_std <<"\n";
    gsInfo << "Difference: \n" << (mat.toSparseMatrix().toDense() - mat_std).cwiseAbs().maxCoeff() <<"\n";

    // rank-1 density with sinusoidal variation
    std::vector<gsMatrix<real_t>> skeleton(2);

    index_t n1 = basis->component(0).size();
    index_t n2 = basis->component(1).size();

    skeleton[0].resize(n1,1);
    skeleton[1].resize(n2,1);

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

    // --- build rho ---
    for (index_t i = 0; i < n1; ++i)
    {
        skeleton[0](i,0) = 1.0 + 0.8 * std::sin(6.0 * M_PI * gx[i]);
    }
    for (index_t j = 0; j < n2; ++j)
    {
       skeleton[1](j,0) = 1.0 + 0.8 * std::cos(6.0 * M_PI * gy[j]);
    }
    // skeleton[1].col(0).setOnes();

    gsTensorFunction<> rho(*basis, skeleton);

    gsStopwatch timer_diffusion;
    timer_diffusion.restart();
    gsTensorAssembler<> ta_diffusion;
    ta_diffusion.compute(*basis, DIFFUSION, rho);
    real_t assembly_time_tensor = timer_diffusion.stop();
    gsInfo << "Diffusion assembly time: " << assembly_time_tensor << " seconds\n";
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

    auto applyDirichletDense = [&](gsMatrix<real_t>& A,
                               gsMatrix<real_t>& rhs,
                               const std::vector<index_t>& dofs,
                               const std::vector<real_t>& vals)
    {
        GISMO_ASSERT(dofs.size() == vals.size(), "Dirichlet size mismatch.");

        for (size_t k = 0; k < dofs.size(); ++k)
        {
            index_t c = dofs[k];
            real_t g  = vals[k];

            // rhs <- rhs - A(:,c) * g
            rhs -= A.col(c) * g;

            // zero row and column
            A.row(c).setZero();
            A.col(c).setZero();

            A(c,c) = 1.0;
            rhs(c,0) = g;
        }
    };

    gsMatrix<real_t> Kx_dense = gsMatrix<real_t>(K);
    gsMatrix<real_t> Ky_dense = gsMatrix<real_t>(K);

    applyDirichletDense(Kx_dense, rhs_x, bdofs, vals_x);
    applyDirichletDense(Ky_dense, rhs_y, bdofs, vals_y);

    gsMatrix<real_t> sol_x = Kx_dense.fullPivLu().solve(rhs_x);
    gsMatrix<real_t> sol_y = Ky_dense.fullPivLu().solve(rhs_y);

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

