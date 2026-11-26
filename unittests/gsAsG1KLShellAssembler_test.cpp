/** @file gsAsG1KLShellAssembler_test.cpp

    @brief Unit tests for gsAsG1KLShellAssembler and gsAsG1ShellHierarchy
           (intrinsic symmetry/rigid-body kernel, physical plate/cantilever/shell
           benchmarks, curved sharp-corner BCs, and validation guards).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "gismo_unittest.h"
#include <gsModeling/gsAsG1ShellSolver.hpp>

using namespace gismo;

namespace {

gsMultiPatch<real_t> embed2DTo3D(const gsMultiPatch<real_t>& mp2)
{
    gsMultiPatch<real_t> mp3;
    for (size_t p = 0; p < mp2.nPatches(); ++p)
    {
        const auto& bsp = dynamic_cast<const gsTensorBSpline<2, real_t>&>(mp2.patch(p));
        gsMatrix<real_t> c3 = gsMatrix<real_t>::Zero(bsp.coefs().rows(), 3);
        c3.leftCols(2) = bsp.coefs();
        mp3.addPatch(typename gsGeometry<real_t>::uPtr(new gsTensorBSpline<2, real_t>(bsp.basis(), c3)));
    }
    mp3.computeTopology();
    return mp3;
}

gsMultiPatch<real_t> makeCurvedParaboloid2x2()
{
    gsMultiPatch<real_t> mp = embed2DTo3D(gsNurbsCreator<real_t>::BSplineSquareGrid(2, 2, 0.5));
    for (size_t p = 0; p < mp.nPatches(); ++p)
    {
        auto& bsp = dynamic_cast<gsTensorBSpline<2, real_t>&>(mp.patch(p));
        bsp.degreeElevate(1);
        const gsMatrix<real_t> anc = bsp.basis().anchors();
        const gsMatrix<real_t> X = bsp.eval(anc);
        gsMatrix<real_t> vals(3, anc.cols());
        for (index_t j = 0; j < anc.cols(); ++j)
        {
            vals(0, j) = X(0, j);
            vals(1, j) = X(1, j);
            vals(2, j) = 0.25 * (X(0, j) * X(0, j) + X(1, j) * X(1, j));
        }
        bsp.coefs() = bsp.basis().interpolateAtAnchors(vals)->coefs();
    }
    mp.computeTopology();
    return mp;
}

} // anonymous namespace

SUITE(gsAsG1KLShellAssembler_test)
{
    TEST(SymmetryRigidBodyAndManufacturedConsistency)
    {
        gsMultiPatch<real_t> mp = makeCurvedParaboloid2x2();
        gsKLShellMaterial<real_t> mat(1e6, 0.3, 0.1);
        gsMultiBasis<real_t> db = buildAsG1PolynomialMultiBasis(mp, 3, 1, 0);
        gsAsG1KLShellAssembler<real_t> asmbl(mp, db, mat);

        gsSparseMatrix<real_t> K;
        asmbl.assembleStiffness(K);
        const real_t kNorm = K.norm();
        CHECK(kNorm > 0.0);

        const real_t symErr = (gsSparseMatrix<real_t>(K - gsSparseMatrix<real_t>(K.transpose()))).norm() / kNorm;
        CHECK(symErr < 1e-14);

        const index_t N = asmbl.scalarDiscontSize();
        for (int m = 0; m < 6; ++m)
        {
            gsVector<real_t> c(3 * N);
            for (size_t p = 0; p < mp.nPatches(); ++p)
            {
                const gsMatrix<real_t> anc = db[p].anchors();
                const gsMatrix<real_t> X = mp.patch(p).eval(anc);
                gsMatrix<real_t> vals(3, anc.cols());
                for (index_t j = 0; j < anc.cols(); ++j)
                {
                    gsVector<real_t, 3> u = gsVector<real_t, 3>::Zero();
                    if (m < 3)
                        u(m) = 1.0;
                    else
                    {
                        gsVector<real_t, 3> ax = gsVector<real_t, 3>::Zero();
                        ax(m - 3) = 1.0;
                        gsVector<real_t, 3> x3;
                        x3 << X(0, j), X(1, j), X(2, j);
                        u = ax.cross(x3);
                    }
                    vals.col(j) = u;
                }
                auto g = db[p].interpolateAtAnchors(vals);
                for (index_t comp = 0; comp < 3; ++comp)
                    c.segment(comp * N + asmbl.patchOffsets()[p], db[p].size()) = g->coefs().col(comp);
            }
            CHECK((K * c).norm() / (kNorm * c.norm()) < 1e-13);
        }

        // Affine displacement field is exactly in the spline space -> F_mfg == K * c
        gsMatrix<real_t, 3, 3> M;
        M << 0.3, -0.2, 0.5, 0.1, 0.4, -0.3, -0.6, 0.2, 0.25;
        gsVector<real_t, 3> b;
        b << 0.1, -0.2, 0.3;
        auto exactAffine = [&](index_t comp, const gsVector<real_t>& X,
                               real_t& v, gsVector<real_t>& g, gsMatrix<real_t>& H) {
            g.setZero(3);
            H.setZero(3, 3);
            v = b(comp);
            for (int j = 0; j < 3; ++j)
            {
                v += M(comp, j) * X(j);
                g(j) = M(comp, j);
            }
        };
        gsVector<real_t> F;
        asmbl.assembleManufacturedRHS(exactAffine, F);
        gsVector<real_t> ca(3 * N);
        for (size_t p = 0; p < mp.nPatches(); ++p)
        {
            const gsMatrix<real_t> anc = db[p].anchors();
            const gsMatrix<real_t> X = mp.patch(p).eval(anc);
            gsMatrix<real_t> vals(3, anc.cols());
            for (index_t j = 0; j < anc.cols(); ++j)
            {
                gsVector<real_t, 3> x3;
                x3 << X(0, j), X(1, j), X(2, j);
                vals.col(j) = M * x3 + b;
            }
            auto g = db[p].interpolateAtAnchors(vals);
            for (index_t comp = 0; comp < 3; ++comp)
                ca.segment(comp * N + asmbl.patchOffsets()[p], db[p].size()) = g->coefs().col(comp);
        }
        CHECK((F - K * ca).norm() / F.norm() < 1e-12);
    }

    TEST(MultiPatchClampedAndSimplySupportedPlate)
    {
        gsMultiPatch<real_t> grid = embed2DTo3D(gsNurbsCreator<real_t>::BSplineSquareGrid(2, 2, 0.5));
        gsKLShellMaterial<real_t> mat(1e6, 0.3, 0.01);
        const real_t D = mat.E * std::pow(mat.thickness, 3) / (12.0 * (1.0 - mat.nu * mat.nu));

        const std::pair<std::string, real_t> cases[2] = {
            {"ValuesAndDerivatives", 0.00126532 / D},
            {"Values",               0.00406235 / D}
        };

        for (const auto& cs : cases)
        {
            gsBoundaryConditions<real_t> bc;
            gsConstantFunction<real_t> zero(0.0, 3);
            for (const auto& ps : grid.boundaries())
                for (short_t c = 0; c < 3; ++c)
                    bc.add(ps.patch, ps.side(), cs.first, zero, c);

            gsAsG1ShellHierarchy<real_t> hier(grid, bc, 3, 1, 0, 1, true);
            const auto& fin = hier.finest();
            gsAsG1KLShellAssembler<real_t> asmbl(grid, fin.dbasis, mat);
            gsSparseMatrix<real_t> K, Kf;
            gsVector<real_t> F, Ff, ub = gsVector<real_t>::Zero(fin.mapper.boundarySize()), uf;
            asmbl.assembleStiffness(K);
            asmbl.assemblePhysicalRHS(
                [](const gsVector<real_t>&, gsVector<real_t, 3>& f, real_t& pn) { f.setZero(); pn = 1.0; }, F);
            asmbl.projectSystem(K, F, fin.T_free, fin.T_bnd, ub, Kf, Ff);
            const auto st = hier.solve(Kf, Ff, mat, uf, "cg", "line_ring1", 1, 1, 1, 200, 1e-12);
            CHECK(st.converged);

            gsVector<real_t> c_sol = fin.T_free * uf;
            real_t wc = 0.0;
            for (size_t p = 0; p < grid.nPatches(); ++p)
                for (int k = 0; k < 4; ++k)
                {
                    gsMatrix<real_t> u(2, 1);
                    u << (k % 2), (k / 2);
                    const gsMatrix<real_t> X = grid.patch(p).eval(u);
                    if (std::abs(X(0, 0) - 0.5) < 1e-10 && std::abs(X(1, 0) - 0.5) < 1e-10)
                        wc = asmbl.evalDisplacement(c_sol, p, u(0, 0), u(1, 0))(2);
                }
            CHECK_CLOSE(1.0, wc / cs.second, 1e-4);
        }
    }

    TEST(MultiPatchCantileverBothOrientations)
    {
        gsMultiPatch<real_t> grid = embed2DTo3D(gsNurbsCreator<real_t>::BSplineSquareGrid(2, 2, 0.5));
        gsKLShellMaterial<real_t> mat(1e6, 0.0, 0.01);
        const real_t D = mat.E * std::pow(mat.thickness, 3) / 12.0;
        const real_t wTipRef = 1.0 / (8.0 * D);

        for (int clampAxis = 0; clampAxis < 2; ++clampAxis)
        {
            gsBoundaryConditions<real_t> bc;
            gsConstantFunction<real_t> zero(0.0, 3);
            for (const auto& ps : grid.boundaries())
            {
                gsMatrix<real_t> pt(2, 1);
                const short_t nDir = ps.direction(), tDir = 1 - nDir;
                pt(nDir, 0) = ps.parameter() ? 1.0 : 0.0;
                pt(tDir, 0) = 0.5;
                const gsMatrix<real_t> X = grid.patch(ps.patch).eval(pt);
                if (std::abs(X(clampAxis, 0)) < 1e-10)
                    for (short_t c = 0; c < 3; ++c)
                        bc.add(ps.patch, ps.side(), "ValuesAndDerivatives", zero, c);
            }

            gsAsG1ShellHierarchy<real_t> hier(grid, bc, 4, 2, 0, 1, true);
            const auto& fin = hier.finest();
            gsAsG1KLShellAssembler<real_t> asmbl(grid, fin.dbasis, mat);
            gsSparseMatrix<real_t> K, Kf;
            gsVector<real_t> F, Ff, ub = gsVector<real_t>::Zero(fin.mapper.boundarySize()), uf;
            asmbl.assembleStiffness(K);
            asmbl.assemblePhysicalRHS(
                [](const gsVector<real_t>&, gsVector<real_t, 3>& f, real_t& pn) { f.setZero(); pn = 1.0; }, F);
            asmbl.projectSystem(K, F, fin.T_free, fin.T_bnd, ub, Kf, Ff);
            const auto st = hier.solve(Kf, Ff, mat, uf, "direct");
            CHECK(st.converged);

            gsVector<real_t> c_sol = fin.T_free * uf;
            for (size_t p = 0; p < grid.nPatches(); ++p)
            {
                for (int j = 0; j <= 4; ++j)
                    for (int i = 0; i <= 4; ++i)
                    {
                        gsMatrix<real_t> pt(2, 1);
                        pt << i * 0.25, j * 0.25;
                        const gsMatrix<real_t> X = grid.patch(p).eval(pt);
                        const auto uVal = asmbl.evalDisplacement(c_sol, p, pt(0, 0), pt(1, 0));
                        if (std::abs(X(clampAxis, 0) - 1.0) < 1e-10)
                            CHECK(std::abs(uVal(2) - wTipRef) / wTipRef < 1e-8);
                        if (std::abs(X(clampAxis, 0)) < 1e-10)
                            CHECK(uVal.norm() < 1e-12);
                    }
            }
        }
    }

    TEST(CurvedSharpCornerOneSidedDirichlet)
    {
        gsMultiPatch<real_t> ann;
        ann.addPatch(gsNurbsCreator<real_t>::NurbsQuarterAnnulus(1.0, 2.0));
        const auto& r2 = dynamic_cast<const gsTensorNurbs<2, real_t>&>(ann.patch(0));
        gsMatrix<real_t> c3 = gsMatrix<real_t>::Zero(r2.coefs().rows(), 3);
        c3.leftCols(2) = r2.coefs();
        gsMultiPatch<real_t> ann3;
        ann3.addPatch(typename gsGeometry<real_t>::uPtr(new gsTensorNurbs<2, real_t>(r2.basis(), c3)));
        ann3.computeTopology();

        gsBoundaryConditions<real_t> bc;
        gsConstantFunction<real_t> zero(0.0, 3);
        std::vector<patchSide> curvedSides;
        for (const auto& ps : ann3.boundaries())
        {
            gsMatrix<real_t> pts(2, 3);
            const short_t nDir = ps.direction(), tDir = 1 - nDir;
            pts.row(nDir).setConstant(ps.parameter() ? 1.0 : 0.0);
            pts(tDir, 0) = 0.0; pts(tDir, 1) = 0.5; pts(tDir, 2) = 1.0;
            const gsMatrix<real_t> X = ann3.patch(0).eval(pts);
            const real_t r0 = std::hypot(X(0, 0), X(1, 0));
            const real_t r1 = std::hypot(X(0, 1), X(1, 1));
            const real_t r2v = std::hypot(X(0, 2), X(1, 2));
            if (std::abs(r0 - r1) < 1e-4 && std::abs(r1 - r2v) < 1e-4)
            {
                curvedSides.push_back(ps);
                for (short_t c = 0; c < 3; ++c)
                    bc.add(ps.patch, ps.side(), "Values", zero, c);
            }
        }
        CHECK_EQUAL(2u, curvedSides.size());

        gsKLShellMaterial<real_t> mat(1e6, 0.3, 0.05);
        gsAsG1ShellHierarchy<real_t> hier(ann3, bc, 4, 2, 0, 1, true);
        const auto& fin = hier.finest();
        gsAsG1KLShellAssembler<real_t> asmbl(ann3, fin.dbasis, mat);
        gsSparseMatrix<real_t> K, Kf;
        gsVector<real_t> F, Ff, ub = gsVector<real_t>::Zero(fin.mapper.boundarySize()), uf;
        asmbl.assembleStiffness(K);
        asmbl.assemblePhysicalRHS(
            [](const gsVector<real_t>&, gsVector<real_t, 3>& f, real_t& pn) { f.setZero(); pn = 1.0; }, F);
        asmbl.projectSystem(K, F, fin.T_free, fin.T_bnd, ub, Kf, Ff);
        const auto st = hier.solve(Kf, Ff, mat, uf, "direct");
        CHECK(st.converged);

        gsVector<real_t> c_sol = fin.T_free * uf;
        for (const patchSide& ps : curvedSides)
        {
            const short_t nDir = ps.direction(), tDir = 1 - nDir;
            for (int s = 0; s <= 20; ++s)
            {
                gsMatrix<real_t> pt(2, 1);
                pt(nDir, 0) = ps.parameter() ? 1.0 : 0.0;
                pt(tDir, 0) = s / 20.0;
                CHECK(asmbl.evalDisplacement(c_sol, 0, pt(0, 0), pt(1, 0)).norm() < 1e-11);
            }
        }
    }

    TEST(ValidationGuards)
    {
        gsMultiPatch<real_t> sq;
        sq.addPatch(gsNurbsCreator<real_t>::BSplineSquare(1.0, 0.0, 0.0));
        sq = embed2DTo3D(sq);

        gsBoundaryConditions<real_t> bc;
        gsConstantFunction<real_t> zero(0.0, 3);
        for (const auto& ps : sq.boundaries())
            for (short_t c = 0; c < 3; ++c)
                bc.add(ps.patch, ps.side(), "ValuesAndDerivatives", zero, c);

        gsKLShellMaterial<real_t> mat(1e6, 0.3, 0.1);
        gsAsG1ShellHierarchy<real_t> hier(sq, bc, 3, 1, 0, 0, true);
        const auto& fin = hier.finest();
        gsAsG1KLShellAssembler<real_t> asmbl(sq, fin.dbasis, mat);

        gsSparseMatrix<real_t> K, Kf;
        gsVector<real_t> F, Ff, uf;
        asmbl.assembleStiffness(K);
        asmbl.assemblePhysicalRHS(
            [](const gsVector<real_t>&, gsVector<real_t, 3>& f, real_t& pn) { f.setZero(); pn = 1.0; }, F);

        // Mis-sized non-empty u_bnd must throw (not silently drop boundary lifting)
        gsVector<real_t> badUb = gsVector<real_t>::Ones(fin.mapper.boundarySize() + 1);
        CHECK_THROW(asmbl.projectSystem(K, F, fin.T_free, fin.T_bnd, badUb, Kf, Ff), std::runtime_error);

        // Empty u_bnd is valid homogeneous shorthand
        gsVector<real_t> emptyUb;
        asmbl.projectSystem(K, F, fin.T_free, fin.T_bnd, emptyUb, Kf, Ff);

        // Unknown solverType must throw
        CHECK_THROW(hier.solve(Kf, Ff, mat, uf, "CG"), std::runtime_error);

        // Mis-sized displacement coefficient vector must throw
        gsVector<real_t> badC = gsVector<real_t>::Zero(10);
        CHECK_THROW(asmbl.evalDisplacement(badC, 0, 0.5, 0.5), std::runtime_error);
    }

    TEST(PhysicalClampedNitscheCurvedShellAndNormalRotation)
    {
        // 1. Curved 2x2 paraboloid physically Clamped shell with nonzero membrane boundary strain
        gsMultiPatch<real_t> mp = makeCurvedParaboloid2x2();
        const real_t alpha = 0.25;
        auto exactClampedCurved = [&](index_t comp, const gsVector<real_t>& X,
                                      real_t& val, gsVector<real_t>& g, gsMatrix<real_t>& H)
        {
            const real_t pi = EIGEN_PI;
            const real_t x = X(0), y = X(1);
            const real_t sx = std::sin(pi * x), cx = std::cos(pi * x);
            const real_t sy = std::sin(pi * y), cy = std::cos(pi * y);
            const real_t phi = sx * sy;
            const real_t px = pi * cx * sy, py = pi * sx * cy;
            const real_t pxx = -pi * pi * phi, pyy = -pi * pi * phi, pxy = pi * pi * cx * cy;

            g.setZero(3);
            H.setZero(3, 3);
            if (comp == 0 || comp == 1)
            {
                val = phi;
                g(0) = px; g(1) = py;
                H(0, 0) = pxx; H(1, 1) = pyy; H(0, 1) = H(1, 0) = pxy;
            }
            else
            {
                const real_t c2 = 2.0 * alpha;
                const real_t gx = sx * sx, gy = sy * sy;
                const real_t dgx = pi * std::sin(2.0 * pi * x), dgy = pi * std::sin(2.0 * pi * y);
                const real_t ddgx = 2.0 * pi * pi * std::cos(2.0 * pi * x);
                const real_t ddgy = 2.0 * pi * pi * std::cos(2.0 * pi * y);
                const real_t psi = gx * gy;
                const real_t wx = dgx * gy, wy = gx * dgy;
                const real_t wxx = ddgx * gy, wyy = gx * ddgy, wxy = dgx * dgy;

                val = c2 * (x + y) * phi + psi;
                g(0) = c2 * phi + c2 * (x + y) * px + wx;
                g(1) = c2 * phi + c2 * (x + y) * py + wy;
                H(0, 0) = 2.0 * c2 * px + c2 * (x + y) * pxx + wxx;
                H(1, 1) = 2.0 * c2 * py + c2 * (x + y) * pyy + wyy;
                H(0, 1) = H(1, 0) = c2 * (px + py) + c2 * (x + y) * pxy + wxy;
            }
        };

        gsBoundaryConditions<real_t> bc;
        gsConstantFunction<real_t> zero(0.0, 3);
        for (const auto& ps : mp.boundaries())
            for (short_t c = 0; c < 3; ++c)
                bc.add(ps.patch, ps.side(), "Clamped", zero, c);

        gsKLShellMaterial<real_t> mat(1e6, 0.3, 0.05);
        real_t l2_0 = 0, h2_0 = 0, rot_0 = 0, l2_1 = 0, h2_1 = 0, rot_1 = 0;
        for (index_t r = 0; r <= 1; ++r)
        {
            gsAsG1ShellHierarchy<real_t> hier(mp, bc, 3, 1, 0, r, true);
            const auto& fin = hier.finest();
            gsAsG1KLShellAssembler<real_t> asmbl(mp, fin.dbasis, mat, bc);
            CHECK_EQUAL(static_cast<index_t>(mp.boundaries().size()), asmbl.numRotationBoundarySides());
            gsSparseMatrix<real_t> K, Kf;
            gsVector<real_t> F, Ff, ub, uf;
            asmbl.assembleStiffness(K);
            asmbl.assembleManufacturedRHS(exactClampedCurved, F);
            asmbl.projectSystem(K, F, fin.T_free, fin.T_bnd, ub, Kf, Ff);
            const auto st = hier.solve(Kf, Ff, mat, uf, "cg", "line_ring1", 1, 1, 1, 200, 1e-12);
            CHECK(st.converged);
            gsVector<real_t> c_sol = fin.T_free * uf;
            const auto err = asmbl.computeErrors(c_sol, exactClampedCurved);
            if (r == 0) { l2_0 = err.l2Rel(); h2_0 = err.h2Rel(); rot_0 = asmbl.computeRotationBCErrorL2(c_sol); }
            else        { l2_1 = err.l2Rel(); h2_1 = err.h2Rel(); rot_1 = asmbl.computeRotationBCErrorL2(c_sol); }
        }
        CHECK(std::log2(l2_0 / l2_1) > 3.5);
        CHECK(std::log2(h2_0 / h2_1) > 1.7);
        CHECK(std::log2(rot_0 / rot_1) > 2.5);
        CHECK(rot_1 < 1e-3);

        // 2. Inhomogeneous NormalRotation patch test (exact tilted plane u_z = slope * x)
        gsMultiPatch<real_t> sq;
        sq.addPatch(gsNurbsCreator<real_t>::BSplineSquare(1.0, 0.0, 0.0));
        sq = embed2DTo3D(sq);
        const real_t slope = 0.05;
        gsBoundaryConditions<real_t> bcRot;
        gsConstantFunction<real_t> rotVal(-slope, 3);
        for (short_t c = 0; c < 3; ++c)
            bcRot.add(0, boundary::west, "Values", zero, c);
        bcRot.add(0, boundary::west, "NormalRotation", rotVal, 0);

        gsAsG1ShellHierarchy<real_t> hSq(sq, bcRot, 3, 1, 0, 0, true);
        const auto& fSq = hSq.finest();
        gsAsG1KLShellAssembler<real_t> aSq(sq, fSq.dbasis, mat, bcRot);
        gsSparseMatrix<real_t> Ksq, Kfsq;
        gsVector<real_t> Fsq, Ffsq, ubsq, ufsq;
        aSq.assembleStiffness(Ksq);
        aSq.assemblePhysicalRHS(
            [](const gsVector<real_t>&, gsVector<real_t, 3>& f, real_t& pn) { f.setZero(); pn = 0.0; }, Fsq);
        aSq.projectSystem(Ksq, Fsq, fSq.T_free, fSq.T_bnd, ubsq, Kfsq, Ffsq);
        const auto stSq = hSq.solve(Kfsq, Ffsq, mat, ufsq, "direct");
        CHECK(stSq.converged);
        gsVector<real_t> cSq = fSq.T_free * ufsq;
        for (int i = 0; i <= 4; ++i)
        {
            const real_t x = i * 0.25;
            const auto u = aSq.evalDisplacement(cSq, 0, x, 0.5);
            CHECK_CLOSE(0.0, u(0), 1e-12);
            CHECK_CLOSE(0.0, u(1), 1e-12);
            CHECK_CLOSE(slope * x, u(2), 1e-12);
        }
        CHECK(aSq.computeRotationBCErrorL2(cSq) < 1e-12);
    }

    TEST(AuditCorrectnessTrapsRegression)
    {
        // F1: verifyBCs catches forgetting bc when "Clamped" is present; makeAssembler binds it
        {
            gsMultiPatch<real_t> mp = makeCurvedParaboloid2x2();
            gsBoundaryConditions<real_t> bc;
            gsConstantFunction<real_t> zero(0.0, 3);
            for (const auto& ps : mp.boundaries())
                for (short_t c = 0; c < 3; ++c)
                    bc.add(ps.patch, ps.side(), "Clamped", zero, c);
            gsAsG1ShellHierarchy<real_t> hier(mp, bc, 3, 1, 0, 0, true);
            gsKLShellMaterial<real_t> mat(1e6, 0.3, 0.05);
            gsAsG1KLShellAssembler<real_t> missingBC(mp, hier.finest().dbasis, mat);
            CHECK_THROW(hier.verifyAssemblerBCs(missingBC), std::runtime_error);

            gsAsG1KLShellAssembler<real_t> boundAsm = hier.makeAssembler(mat);
            hier.verifyAssemblerBCs(boundAsm);

            // Unsupported label rejected loudly
            gsBoundaryConditions<real_t> badBC = bc;
            badBC.add(mp.boundaries()[0].patch, mp.boundaries()[0].side(), "Neumann", zero, 2);
            CHECK_THROW(boundAsm.setRotationBCs(badBC), std::runtime_error);
        }

        // F2: Inhomogeneous physical Dirichlet lifting reproduces exact rigid translation on curved clamped shell
        {
            gsMultiPatch<real_t> mp = makeCurvedParaboloid2x2();
            gsVector<real_t> shift(3);
            shift << 1.5, -0.7, 2.3;
            gsConstantFunction<real_t> shiftFn(shift, 3);
            gsBoundaryConditions<real_t> bc;
            for (const auto& ps : mp.boundaries())
                for (short_t c = 0; c < 3; ++c)
                    bc.add(ps.patch, ps.side(), "Clamped", shiftFn, c);

            gsKLShellMaterial<real_t> mat(1e6, 0.5, 0.05); // Also verifies F13a (nu == 0.5 incompressible plane stress)
            gsAsG1ShellHierarchy<real_t> hier(mp, bc, 3, 1, 0, 1, true);
            gsAsG1KLShellAssembler<real_t> asmbl = hier.makeAssembler(mat);
            gsSparseMatrix<real_t> K, Kf;
            gsVector<real_t> F, Ff, ub, uf;
            asmbl.assembleStiffness(K);
            asmbl.assemblePhysicalRHS(
                [](const gsVector<real_t>&, gsVector<real_t, 3>& f, real_t& pn) { f.setZero(); pn = 0.0; }, F);
            hier.projectSystem(asmbl, K, F, Kf, Ff, &ub);
            CHECK(ub.size() == hier.finest().mapper.boundarySize());
            CHECK(ub.norm() > 0.0);

            const auto st = hier.solve(Kf, Ff, mat, uf, "cg", "line_ring1", 1, 1, 1, 200, 1e-12);
            CHECK(st.converged);
            gsVector<real_t> c_sol = hier.finest().T_free * uf + hier.finest().T_bnd * ub;
            for (size_t p = 0; p < mp.nPatches(); ++p)
            {
                const auto uMid = asmbl.evalDisplacement(c_sol, p, 0.5, 0.5);
                CHECK_CLOSE(shift(0), uMid(0), 1e-11);
                CHECK_CLOSE(shift(1), uMid(1), 1e-11);
                CHECK_CLOSE(shift(2), uMid(2), 1e-11);
            }
            CHECK(asmbl.computeRotationBCErrorL2(c_sol) < 1e-11);
        }

        // F3: Creased 45-degree 2-patch folded L-shell must be rejected by computeGluingData
        {
            gsMultiPatch<real_t> folded;
            gsTensorBSpline<2, real_t> p0 = *dynamic_cast<gsTensorBSpline<2, real_t>*>(
                embed2DTo3D(gsMultiPatch<real_t>(*gsNurbsCreator<real_t>::BSplineSquare(1.0, 0.0, 0.0))).patch(0).clone().release());
            gsTensorBSpline<2, real_t> p1 = p0;
            const real_t c45 = std::cos(EIGEN_PI / 4.0), s45 = std::sin(EIGEN_PI / 4.0);
            for (index_t i = 0; i < p1.coefs().rows(); ++i)
            {
                const real_t x = p1.coefs()(i, 0);
                p1.coefs()(i, 0) = 1.0 + x * c45;
                p1.coefs()(i, 2) = x * s45;
            }
            folded.addPatch(p0);
            folded.addPatch(p1);
            folded.computeTopology();
            CHECK_THROW(computeGluingData(folded, real_t(1e-5)), std::runtime_error);
        }

        // F5 & F6: Acutely skewed 20-deg rhombus clamped plate: SymmetricNitsche is SPD & accurate,
        //          while componentwise "Derivatives" on a non-orthogonal sharp corner throws cleanly.
        {
            const real_t ang = 20.0 * EIGEN_PI / 180.0;
            gsMultiPatch<real_t> rhomb = embed2DTo3D(gsMultiPatch<real_t>(*gsNurbsCreator<real_t>::BSplineSquare(1.0, 0.0, 0.0)));
            auto& bsp = dynamic_cast<gsTensorBSpline<2, real_t>&>(rhomb.patch(0));
            for (index_t i = 0; i < bsp.coefs().rows(); ++i)
            {
                const real_t u = bsp.coefs()(i, 0), v = bsp.coefs()(i, 1);
                bsp.coefs()(i, 0) = u + v * std::cos(ang);
                bsp.coefs()(i, 1) = v * std::sin(ang);
            }
            rhomb.computeTopology();

            gsConstantFunction<real_t> zero(0.0, 3);
            // F6: "Derivatives" on one side of non-orthogonal sharp corner must throw
            gsBoundaryConditions<real_t> bcBadDer;
            bcBadDer.add(0, boundary::west, "Derivatives", zero, 2);
            CHECK_THROW(gsAsG1ShellHierarchy<real_t>(rhomb, bcBadDer, 3, 1, 0, 0, true), std::runtime_error);

            // F5: Clamped 20-deg rhombus with element-local eigenvalue SymmetricNitsche is strictly SPD
            gsBoundaryConditions<real_t> bcClamp;
            for (const auto& ps : rhomb.boundaries())
                for (short_t c = 0; c < 3; ++c)
                    bcClamp.add(ps.patch, ps.side(), "Clamped", zero, c);

            gsKLShellMaterial<real_t> mat(1e6, 0.3, 0.05);
            gsAsG1ShellHierarchy<real_t> hier(rhomb, bcClamp, 3, 1, 0, 1, true);
            gsAsG1KLShellAssembler<real_t> asmbl = hier.makeAssembler(mat);
            gsSparseMatrix<real_t> K, Kf;
            gsVector<real_t> F, Ff, uf;
            asmbl.assembleStiffness(K);
            asmbl.assemblePhysicalRHS(
                [](const gsVector<real_t>&, gsVector<real_t, 3>& f, real_t& pn) { f.setZero(); pn = 1.0; }, F);
            hier.projectSystem(asmbl, K, F, Kf, Ff);

            gsEigen::SimplicialLDLT<gsSparseMatrix<real_t>> ldlt;
            ldlt.compute(Kf);
            CHECK(ldlt.info() == gsEigen::Success);
            CHECK(ldlt.vectorD().minCoeff() > 0.0);

            const auto st = hier.solve(Kf, Ff, mat, uf, "cg", "line_ring1", 1, 1, 1, 200, 1e-11);
            CHECK(st.converged);
            gsVector<real_t> c_sol = hier.finest().T_free * uf;
            CHECK(asmbl.evalDisplacement(c_sol, 0, 0.5, 0.5)(2) > 0.0);

            // F10: flipGlobalNormalOrientation negates pressure-driven deflection exactly
            asmbl.flipGlobalNormalOrientation();
            gsVector<real_t> Fflip, Ffflip, ufflip;
            asmbl.assemblePhysicalRHS(
                [](const gsVector<real_t>&, gsVector<real_t, 3>& f, real_t& pn) { f.setZero(); pn = 1.0; }, Fflip);
            hier.projectSystem(asmbl, K, Fflip, Kf, Ffflip);
            hier.solve(Kf, Ffflip, mat, ufflip, "direct");
            CHECK_CLOSE(-uf.norm(), -ufflip.norm(), 1e-12);
            CHECK((uf + ufflip).norm() <= 1e-11 * uf.norm());

            // F12: collectShellCoupledStripDofs with r_p = -1 activates default p + 1 rings
            const auto skelNeg = shell_internal::collectShellCoupledStripDofs(
                rhomb, hier.finest().dbasis, hier.finest().argBasis, hier.finest().mapper, -1, true, 3);
            const auto skelZero = shell_internal::collectShellCoupledStripDofs(
                rhomb, hier.finest().dbasis, hier.finest().argBasis, hier.finest().mapper, 0, true, 3);
            CHECK(skelNeg.size() > skelZero.size());
        }

        // F8: Mismatched tangential parametric lengths across interface must throw
        {
            gsMultiPatch<real_t> badParam = embed2DTo3D(gsNurbsCreator<real_t>::BSplineSquareGrid(2, 1, 0.5));
            auto& bsp1 = dynamic_cast<gsTensorBSpline<2, real_t>&>(badParam.patch(1));
            bsp1.basis().knots(1).affineTransformTo(0.0, 2.0);
            badParam.computeTopology();
            CHECK_THROW(computeGluingData(badParam, real_t(1e-5)), std::runtime_error);
        }

        // F9: Missing interior geometry breakpoint in dbasis must throw
        {
            gsMultiPatch<real_t> mp = embed2DTo3D(gsMultiPatch<real_t>(*gsNurbsCreator<real_t>::BSplineSquare(1.0, 0.0, 0.0)));
            mp.patch(0).uniformRefine(1); // geometry now has interior knot at 0.5
            gsKnotVector<real_t> kv(0.0, 1.0, 0, 4); // basis missing knot at 0.5
            gsMultiBasis<real_t> badBasis;
            badBasis.addBasis(new gsTensorBSplineBasis<2, real_t>(kv, kv));
            gsKLShellMaterial<real_t> mat(1e6, 0.3, 0.05);
            CHECK_THROW(gsAsG1KLShellAssembler<real_t>(mp, badBasis, mat), std::runtime_error);
        }
    }
}
