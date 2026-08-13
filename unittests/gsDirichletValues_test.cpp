/** @file gsDirichletValues_test.cpp

    @brief Tests gsDirichletValuesByL2Projection (dirichlet::l2Projection)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
*/

#include "gismo_unittest.h"

#include <string>
#include <vector>

namespace
{

// Compare two Dirichlet fixed-DoF vectors keyed by (patch, local dof, component)
// instead of by boundary index, so the comparison is invariant under any
// renumbering of the boundary indices by gsDofMapper::eliminateDof.
// Returns max |a - b| over the compared entries; nCompared counts them.
real_t bdofMaxDiff(const gsMultiBasis<real_t> & mb,
                   const std::vector<patchSide> & sides,
                   const std::vector<index_t>   & comps,
                   const gsDofMapper & mapA, const gsMatrix<real_t> & fxA,
                   const gsDofMapper & mapB, const gsMatrix<real_t> & fxB,
                   index_t & nCompared, index_t & nSkipped)
{
    real_t d = 0.0;
    nCompared = 0;
    nSkipped  = 0;
    for (size_t s = 0; s != sides.size(); ++s)
    {
        const index_t pi = sides[s].patch;
        const gsMatrix<index_t> bnd = mb.basis(pi).boundary(sides[s].side());
        for (size_t ci = 0; ci != comps.size(); ++ci)
        {
            const index_t c = comps[ci];
            for (index_t l = 0; l != bnd.size(); ++l)
            {
                const index_t glA = mapA.index(bnd.at(l), pi, c);
                const index_t glB = mapB.index(bnd.at(l), pi, c);
                // MANDATORY guard: bindex/global_to_bindex are documented as
                // undefined for a dof that is not on the boundary.
                // A skip here means a dof you expected to be eliminated was NOT
                // (e.g. a bc silently dropped) -- counted, never swallowed.
                if (!mapA.is_boundary_index(glA) || !mapB.is_boundary_index(glB))
                { ++nSkipped; continue; }
                d = math::max(d, math::abs( fxA.at(mapA.global_to_bindex(glA))
                                          - fxB.at(mapB.global_to_bindex(glB)) ));
                ++nCompared;
            }
        }
    }
    return d;
}

} // anonymous namespace

SUITE(gsDirichletValues_test)
{

    // Dirichlet data lying exactly in the spline space: the L2 projection
    // must reproduce it, hence equal dirichlet::interpolation. Also carries
    // a permanent negative control: L2-projecting data that is one degree
    // too high (out of the space) must differ from the in-space reference
    // by more than a comfortable margin, proving the comparison is live.
    TEST(ReproductionSquare2D)
    {
        gsMultiPatch<> mp;
        mp.addPatch(gsNurbsCreator<>::BSplineSquare(1.0));
        mp.computeTopology();

        const std::vector<patchSide> sides = mp.boundaries();
        const std::vector<index_t> comps(1, 0);

        for (index_t p = 1; p <= 3; ++p)
        {
            gsMultiBasis<> mb(mp);
            if (p > 1) mb.degreeElevate(p-1);
            mb.uniformRefine(3); // 4 elements/direction

            const std::string P = std::to_string(p);
            gsFunctionExpr<> g    ("(1+x)^"+P+"*(1+y)^"+P, 2);
            // Out-of-space perturbation: one degree higher than the basis in x.
            gsFunctionExpr<> gPert("(1+x)^"+P+"*(1+y)^"+P+"+x^"+std::to_string(p+1), 2);

            gsBoundaryConditions<> bc;
            bc.setGeoMap(mp);
            for (auto & bs : mp.boundaries())
                bc.addCondition(bs.patch, bs.side(), condition_type::dirichlet, &g, 0, false, -1);

            gsBoundaryConditions<> bcPert;
            bcPert.setGeoMap(mp);
            for (auto & bs : mp.boundaries())
                bcPert.addCondition(bs.patch, bs.side(), condition_type::dirichlet, &gPert, 0, false, -1);

            gsExprAssembler<> A(1,1);
            A.setIntegrationElements(mb);
            auto u = A.getSpace(mb);

            u.setup(bc, dirichlet::l2Projection, 0);
            const gsMatrix<real_t> fxL2  = u.fixedPart(); // deep copy
            const gsDofMapper      mapL2 = u.mapper();    // deep copy

            u.setup(bc, dirichlet::interpolation, 0);
            const gsMatrix<real_t> fxIn  = u.fixedPart();
            const gsDofMapper      mapIn = u.mapper();

            u.setup(bcPert, dirichlet::l2Projection, 0);
            const gsMatrix<real_t> fxL2Pert  = u.fixedPart();
            const gsDofMapper      mapL2Pert = u.mapper();

            const real_t scale = math::max(fxIn.cwiseAbs().maxCoeff(), (real_t)1);

            index_t nCompared = 0, nSkipped = 0;
            const real_t relErr = bdofMaxDiff(mb, sides, comps,
                                               mapL2, fxL2, mapIn, fxIn,
                                               nCompared, nSkipped) / scale;
            CHECK_CLOSE(relErr, 0.0, 1e-10);

            CHECK(mapL2.boundarySize() > 0);
            CHECK_EQUAL(mapL2.boundarySize(), fxL2.rows());
            CHECK_EQUAL(mapIn.boundarySize(), fxIn.rows());
            CHECK(fxL2.cwiseAbs().maxCoeff() > 1e-8);
            CHECK(nCompared > 0);
            CHECK_EQUAL(0, nSkipped);
            CHECK(fxL2.allFinite());

            // Permanent negative control: perturbed (out-of-space) L2 result
            // must differ from the in-space reference well beyond tolerance.
            index_t nComparedPert = 0, nSkippedPert = 0;
            const real_t relErrPert = bdofMaxDiff(mb, sides, comps,
                                                   mapL2Pert, fxL2Pert, mapIn, fxIn,
                                                   nComparedPert, nSkippedPert) / scale;
            CHECK_EQUAL(0, nSkippedPert);
            CHECK(relErrPert > 1e-6);
        }
    }

    // 3D counterpart of ReproductionSquare2D.
    TEST(ReproductionCube3D)
    {
        gsMultiPatch<> mp;
        mp.addPatch(gsNurbsCreator<>::BSplineCube(1.0));
        mp.computeTopology();

        const std::vector<patchSide> sides = mp.boundaries();
        const std::vector<index_t> comps(1, 0);

        for (index_t p = 1; p <= 3; ++p)
        {
            gsMultiBasis<> mb(mp);
            if (p > 1) mb.degreeElevate(p-1);
            mb.uniformRefine(1); // 2 elements/direction

            const std::string P = std::to_string(p);
            gsFunctionExpr<> g    ("(1+x)^"+P+"*(1+y)^"+P+"*(1+z)^"+P, 3);
            // Out-of-space perturbation: one degree higher than the basis in x.
            gsFunctionExpr<> gPert("(1+x)^"+P+"*(1+y)^"+P+"*(1+z)^"+P+"+x^"+std::to_string(p+1), 3);

            gsBoundaryConditions<> bc;
            bc.setGeoMap(mp);
            for (auto & bs : mp.boundaries())
                bc.addCondition(bs.patch, bs.side(), condition_type::dirichlet, &g, 0, false, -1);

            gsBoundaryConditions<> bcPert;
            bcPert.setGeoMap(mp);
            for (auto & bs : mp.boundaries())
                bcPert.addCondition(bs.patch, bs.side(), condition_type::dirichlet, &gPert, 0, false, -1);

            gsExprAssembler<> A(1,1);
            A.setIntegrationElements(mb);
            auto u = A.getSpace(mb);

            u.setup(bc, dirichlet::l2Projection, 0);
            const gsMatrix<real_t> fxL2  = u.fixedPart();
            const gsDofMapper      mapL2 = u.mapper();

            u.setup(bc, dirichlet::interpolation, 0);
            const gsMatrix<real_t> fxIn  = u.fixedPart();
            const gsDofMapper      mapIn = u.mapper();

            u.setup(bcPert, dirichlet::l2Projection, 0);
            const gsMatrix<real_t> fxL2Pert  = u.fixedPart();
            const gsDofMapper      mapL2Pert = u.mapper();

            const real_t scale = math::max(fxIn.cwiseAbs().maxCoeff(), (real_t)1);

            index_t nCompared = 0, nSkipped = 0;
            const real_t relErr = bdofMaxDiff(mb, sides, comps,
                                               mapL2, fxL2, mapIn, fxIn,
                                               nCompared, nSkipped) / scale;
            CHECK_CLOSE(relErr, 0.0, 1e-10);

            CHECK(mapL2.boundarySize() > 0);
            CHECK_EQUAL(mapL2.boundarySize(), fxL2.rows());
            CHECK_EQUAL(mapIn.boundarySize(), fxIn.rows());
            CHECK(fxL2.cwiseAbs().maxCoeff() > 1e-8);
            CHECK(nCompared > 0);
            CHECK_EQUAL(0, nSkipped);
            CHECK(fxL2.allFinite());

            index_t nComparedPert = 0, nSkippedPert = 0;
            const real_t relErrPert = bdofMaxDiff(mb, sides, comps,
                                                   mapL2Pert, fxL2Pert, mapIn, fxIn,
                                                   nComparedPert, nSkippedPert) / scale;
            CHECK_EQUAL(0, nSkippedPert);
            CHECK(relErrPert > 1e-6);
        }
    }

    // Vector-valued space (u.dim() == 3) on an 8-patch grid: exercises the
    // per-r loop and both the comp != -1 and comp == -1 branches.
    TEST(MultiPatchVectorValued)
    {
        // BSplineCubeGrid already calls computeTopology() internally.
        gsMultiPatch<> mp = gsNurbsCreator<>::BSplineCubeGrid(2,2,2,1.0);

        gsMultiBasis<> mb(mp);
        mb.degreeElevate(1);   // p = 2
        mb.uniformRefine(1);   // one refine, pinned per spec

        // Three distinct scalar data functions of per-direction degree <= 2.
        gsFunctionExpr<> g0("1+x^2*y^2", 3);
        gsFunctionExpr<> g1("2+y^2*z^2", 3);
        gsFunctionExpr<> g2("3+x^2*z^2", 3);
        gsFunctionExpr<> gvec("1+x^2*y^2","2+y^2*z^2","3+x^2*z^2", 3);

        // Out-of-space perturbation of component 0 only (degree 3 in x).
        gsFunctionExpr<> g0Pert("1+x^2*y^2+x^3", 3);

        gsBoundaryConditions<> bcA;
        bcA.setGeoMap(mp);
        for (auto & bs : mp.boundaries())
        {
            // Same (patch,side,unknown) with distinct comp: addCondition dedupes
            // via boundary_condition::isSame, which compares component too, so
            // these are three distinct, all-kept conditions.
            bcA.addCondition(bs.patch, bs.side(), condition_type::dirichlet, &g0, 0, false, 0);
            bcA.addCondition(bs.patch, bs.side(), condition_type::dirichlet, &g1, 0, false, 1);
            bcA.addCondition(bs.patch, bs.side(), condition_type::dirichlet, &g2, 0, false, 2);
        }

        gsBoundaryConditions<> bcAPert;
        bcAPert.setGeoMap(mp);
        for (auto & bs : mp.boundaries())
        {
            bcAPert.addCondition(bs.patch, bs.side(), condition_type::dirichlet, &g0Pert, 0, false, 0);
            bcAPert.addCondition(bs.patch, bs.side(), condition_type::dirichlet, &g1,     0, false, 1);
            bcAPert.addCondition(bs.patch, bs.side(), condition_type::dirichlet, &g2,     0, false, 2);
        }

        gsBoundaryConditions<> bcB;
        bcB.setGeoMap(mp);
        for (auto & bs : mp.boundaries())
            bcB.addCondition(bs.patch, bs.side(), condition_type::dirichlet, &gvec, 0, false, -1);

        gsExprAssembler<> A(1,1);
        A.setIntegrationElements(mb);
        auto u = A.getSpace(mb, 3);

        u.setup(bcA, dirichlet::l2Projection, 0);
        const gsMatrix<real_t> fxA_L2  = u.fixedPart();
        const gsDofMapper      mapA_L2 = u.mapper();

        u.setup(bcA, dirichlet::interpolation, 0);
        const gsMatrix<real_t> fxA_In  = u.fixedPart();
        const gsDofMapper      mapA_In = u.mapper();

        u.setup(bcB, dirichlet::l2Projection, 0);
        const gsMatrix<real_t> fxB_L2  = u.fixedPart();
        const gsDofMapper      mapB_L2 = u.mapper();

        u.setup(bcAPert, dirichlet::l2Projection, 0);
        const gsMatrix<real_t> fxAPert_L2  = u.fixedPart();
        const gsDofMapper      mapAPert_L2 = u.mapper();

        // NOTE: dirichlet::interpolation is NOT a valid oracle for bcB (comp == -1,
        // u.dim() > 1): gsDirichletValuesByTPInterpolation reads dVals.at(l), which
        // is column-major linear indexing into an (h->size() x targetDim) matrix, so
        // for every component r it reads column 0 only (gsDirichletValues.h:157-165).
        // That is a latent library bug (confounder (a) in the task spec), out of
        // scope here -- bcB interpolation is deliberately never computed or compared.

        const std::vector<patchSide> sides = mp.boundaries();
        const std::vector<index_t> compsAll{0,1,2};

        const real_t scale = math::max(fxA_In.cwiseAbs().maxCoeff(), (real_t)1);

        // Check 1: bcA L2 vs bcA interpolation -- valid oracle, comp != -1 everywhere.
        index_t nCompared1 = 0, nSkipped1 = 0;
        const real_t relErr1 = bdofMaxDiff(mb, sides, compsAll,
                                            mapA_L2, fxA_L2, mapA_In, fxA_In,
                                            nCompared1, nSkipped1) / scale;
        CHECK_CLOSE(relErr1, 0.0, 1e-10);

        CHECK(mapA_L2.boundarySize() > 0);
        CHECK_EQUAL(mapA_L2.boundarySize(), fxA_L2.rows());
        CHECK_EQUAL(mapA_In.boundarySize(), fxA_In.rows());
        CHECK(fxA_L2.cwiseAbs().maxCoeff() > 1e-8);
        CHECK(nCompared1 > 0);
        CHECK_EQUAL(0, nSkipped1);
        CHECK(fxA_L2.allFinite());

        // Check 2: bcB L2 (com == -1 / per-r branch) vs bcA L2 (three comp != -1
        // conditions) -- mathematically identical assembly, same matrix per r.
        index_t nCompared2 = 0, nSkipped2 = 0;
        const real_t relErr2 = bdofMaxDiff(mb, sides, compsAll,
                                            mapB_L2, fxB_L2, mapA_L2, fxA_L2,
                                            nCompared2, nSkipped2) / scale;
        CHECK_CLOSE(relErr2, 0.0, 1e-10);

        CHECK(mapB_L2.boundarySize() > 0);
        CHECK_EQUAL(mapB_L2.boundarySize(), fxB_L2.rows());
        CHECK(fxB_L2.cwiseAbs().maxCoeff() > 1e-8);
        CHECK(nCompared2 > 0);
        CHECK_EQUAL(0, nSkipped2);
        CHECK(fxB_L2.allFinite());

        // Permanent negative control: perturbing component 0 out of the space
        // must be visible against the in-space reference.
        index_t nComparedPert = 0, nSkippedPert = 0;
        const real_t relErrPert = bdofMaxDiff(mb, sides, compsAll,
                                               mapAPert_L2, fxAPert_L2, mapA_In, fxA_In,
                                               nComparedPert, nSkippedPert) / scale;
        CHECK_EQUAL(0, nSkippedPert);
        CHECK(relErrPert > 1e-6);
    }

    // Single-patch, two adjacent Dirichlet faces (west, south) sharing an edge:
    // the number of active boundary functions per element (eltBdryFcts.size())
    // varies -- 9 for face-interior elements, 15 for the one element touching
    // both faces (p = 2: 3x3 layer per face, 3 dofs shared on the edge). With
    // fewer than 2 elements/direction every west-face element also touches the
    // south face and n never varies, defeating the point of the test.
    TEST(AdjacentFacesVaryingActives)
    {
        gsMultiPatch<> mp;
        mp.addPatch(gsNurbsCreator<>::BSplineCube(1.0));
        mp.computeTopology();

        gsMultiBasis<> mb(mp);
        mb.degreeElevate(1);   // p = 2
        mb.uniformRefine(3);   // 4 elements/direction

        gsFunctionExpr<> g    ("(1+x)^2*(1+y)^2*(1+z)^2", 3);
        // Out-of-space perturbation (degree 3 in x, basis is degree 2).
        gsFunctionExpr<> gPert("(1+x)^2*(1+y)^2*(1+z)^2+x^3", 3);

        gsBoundaryConditions<> bc;
        bc.setGeoMap(mp);
        bc.addCondition(0, boundary::west,  condition_type::dirichlet, &g, 0, false, -1);
        bc.addCondition(0, boundary::south, condition_type::dirichlet, &g, 0, false, -1);

        gsBoundaryConditions<> bcPert;
        bcPert.setGeoMap(mp);
        bcPert.addCondition(0, boundary::west,  condition_type::dirichlet, &gPert, 0, false, -1);
        bcPert.addCondition(0, boundary::south, condition_type::dirichlet, &gPert, 0, false, -1);

        gsExprAssembler<> A(1,1);
        A.setIntegrationElements(mb);
        auto u = A.getSpace(mb);

        u.setup(bc, dirichlet::l2Projection, 0);
        const gsMatrix<real_t> fxL2  = u.fixedPart();
        const gsDofMapper      mapL2 = u.mapper();

        u.setup(bc, dirichlet::interpolation, 0);
        const gsMatrix<real_t> fxIn  = u.fixedPart();
        const gsDofMapper      mapIn = u.mapper();

        u.setup(bcPert, dirichlet::l2Projection, 0);
        const gsMatrix<real_t> fxL2Pert  = u.fixedPart();
        const gsDofMapper      mapL2Pert = u.mapper();

        std::vector<patchSide> sides;
        sides.push_back(patchSide(0, boundary::west));
        sides.push_back(patchSide(0, boundary::south));
        const std::vector<index_t> comps(1, 0);

        const real_t scale = math::max(fxIn.cwiseAbs().maxCoeff(), (real_t)1);

        // Every eliminated DoF is covered by west or south, so the projection
        // matrix is SPD; a non-finite result here is the finding, not a
        // tolerance issue.
        index_t nCompared = 0, nSkipped = 0;
        const real_t relErr = bdofMaxDiff(mb, sides, comps,
                                           mapL2, fxL2, mapIn, fxIn,
                                           nCompared, nSkipped) / scale;
        CHECK_CLOSE(relErr, 0.0, 1e-10);

        CHECK(mapL2.boundarySize() > 0);
        CHECK_EQUAL(mapL2.boundarySize(), fxL2.rows());
        CHECK_EQUAL(mapIn.boundarySize(), fxIn.rows());
        CHECK(fxL2.cwiseAbs().maxCoeff() > 1e-8);
        CHECK(nCompared > 0);
        CHECK_EQUAL(0, nSkipped);
        CHECK(fxL2.allFinite());

        index_t nComparedPert = 0, nSkippedPert = 0;
        const real_t relErrPert = bdofMaxDiff(mb, sides, comps,
                                               mapL2Pert, fxL2Pert, mapIn, fxIn,
                                               nComparedPert, nSkippedPert) / scale;
        CHECK_EQUAL(0, nSkippedPert);
        CHECK(relErrPert > 1e-6);
    }

}
