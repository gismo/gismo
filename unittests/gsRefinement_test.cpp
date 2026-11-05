/** @file gsRefinement_test.cpp

    @brief Tests refinement

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither
**/

#include "gismo_unittest.h"
#include <random>


using namespace gismo;

bool compareKV(const gsKnotVector<>& kv1, const gsKnotVector<>& kv2)
{
    return kv1.size() == kv2.size()
        && std::equal(kv1.begin(), kv1.end(), kv2.begin());
}


void testBoehm_helper(const gsBSpline<>& bsp, const std::vector<real_t> knots)
{
    const gsKnotVector<> kv_orig = bsp.knots();
    const gsMatrix<> coef_orig = bsp.coefs();

    gsKnotVector<> kv1 = kv_orig;
    gsMatrix<> coef1 = coef_orig;
    for (size_t i = 0; i < knots.size(); ++i)
        gsBoehm(kv1, coef1, knots[i]);

    gsKnotVector<> kv2 = kv_orig;
    gsMatrix<> coef2 = coef_orig;
    gsBoehmRefine(kv2, coef2, kv2.degree(), knots.begin(), knots.end());

    CHECK (compareKV(kv1, kv2));
    CHECK ((coef1 - coef2).array().abs().maxCoeff() <= 1e-12);
}

// copied from gsNorms.hpp, which is no longer in stable
template <typename T>
T computeMaximumDistance(const gsFunction<T>& f1, const gsFunction<T>& f2, const gsVector<T>& lower, const gsVector<T>& upper, int numSamples=1000)
{
    GISMO_ASSERT( f1.domainDim() == f2.domainDim(), "Functions need to have same domain dimension");
    GISMO_ASSERT( f1.targetDim() == f2.targetDim(), "Functions need to have same target dimension");

    gsMatrix<T> points = uniformPointGrid(lower, upper, numSamples);

    gsMatrix<T> values1, values2;
    f1.eval_into( points, values1 );
    f2.eval_into( points, values2 );

    return (values1 - values2).array().abs().maxCoeff();
}

// Round-off check of gsBSpline::degreeElevate: the curve error is measured in
// units of eps * max|c| on 2001 uniform samples of [0,1].
static real_t degreeElevateRelError(const gsKnotVector<> & kv, const gsMatrix<> & coefs,
                                    short_t i, bool & shapeOk)
{
    gsBSpline<> c(kv, coefs);
    gsBSpline<> e(c);
    e.degreeElevate(i);

    gsBSplineBasis<> bref(kv);
    bref.degreeElevate(i);
    shapeOk = ( e.degree() == c.degree() + i ) && ( e.basis().size() == bref.size() );

    const real_t eps = std::numeric_limits<real_t>::epsilon();
    gsVector<> lo(1), up(1);
    lo[0] = 0; up[0] = 1;
    const real_t err = computeMaximumDistance<real_t>(c, e, lo, up, 2001);
    return err / (eps * c.coefs().cwiseAbs().maxCoeff());
}

SUITE(gsRefinement_test)
{
    TEST(testDegreeElevateAccuracy)
    {
        std::mt19937 gen(5489u);
        const real_t bound = 1000;
        for (short_t p = 0; p <= 6; ++p)
            for (short_t i = 1; i <= 2; ++i)
            {
                const index_t sizes[4] = {p + 1, 20, 100, 400};
                for (index_t n : sizes)
                {
                    gsKnotVector<> kv(0, 1, n - p - 1, p + 1);
                    gsMatrix<> coefs(n, 2);
                    for (index_t k = 0; k < coefs.size(); ++k)
                        coefs(k) = 2.0 * gen() / 4294967295.0 - 1.0;
                    bool shapeOk = false;
                    const real_t rel = degreeElevateRelError(kv, coefs, i, shapeOk);
                    CHECK(shapeOk);
                    CHECK(rel <= bound);
                    if (rel > bound)
                        gsInfo << "p=" << p << " i=" << i << " n=" << n << " err/(eps*max|c|)=" << rel << "\n";
                }
            }

        // interior multiplicities cycling 1,..,p
        const short_t p = 6;
        for (short_t i = 1; i <= 2; ++i)
        {
            gsKnotVector<> kv(0, 1, 18, p + 1);
            const std::vector<real_t> bp = kv.unique();
            for (size_t k = 1; k + 1 < bp.size(); ++k)
                kv.insert(bp[k], static_cast<short_t>((k - 1) % p));
            gsMatrix<> coefs(kv.size() - p - 1, 2);
            for (index_t k = 0; k < coefs.size(); ++k)
                coefs(k) = 2.0 * gen() / 4294967295.0 - 1.0;
            bool shapeOk = false;
            const real_t rel = degreeElevateRelError(kv, coefs, i, shapeOk);
            CHECK(shapeOk);
            CHECK(rel <= bound);
            if (rel > bound)
                gsInfo << "multiplicity case i=" << i << " err/(eps*max|c|)=" << rel << "\n";
        }
    }

    TEST(testBoehm)
    {
        UnitTest::deactivate_output();
        gsGeometry<>::uPtr geo( gsNurbsCreator<>::BSplineFatCircle() );
        UnitTest::reactivate_output();
        gsBSpline<>& bsp = dynamic_cast<gsBSpline<>&>(*geo);

        // get unique knots
        std::vector<real_t> knots = bsp.basis().knots().unique();
        // insert a few knots
        const int n = knots.size() - 1;
        for (int i = 0; i < n; ++i)
            for (int j = 1; j <= 3; ++j)
                knots.push_back((4-j)/4.0 * knots[i] + j/4.0 * knots[i+1]);
        std::sort(knots.begin(), knots.end());
        // remove start and end knot
        knots.erase(knots.begin());
        knots.pop_back();

        testBoehm_helper(bsp, knots);
    }

    TEST(testUniformRefine)
    {
        UnitTest::deactivate_output();
        gsGeometry<>::uPtr geo( gsNurbsCreator<>::BSplineFatCircle() );
        UnitTest::reactivate_output();
        gsBSpline<>& bsp = dynamic_cast<gsBSpline<>&>(*geo);
        gsBSpline<>::uPtr bsp_ref = bsp.clone();
        bsp_ref->uniformRefine();

        gsMatrix<> paramRange = bsp.parameterRange();
        const real_t dist = computeMaximumDistance<real_t>(bsp, *bsp_ref, paramRange.col(0), paramRange.col(1));

        CHECK (dist <= 1e-14);
    }

    TEST(testRefineWithMultiplicity)
    {
        UnitTest::deactivate_output();
        gsGeometry<>::uPtr geo( gsNurbsCreator<>::BSplineFatCircle() );
        UnitTest::reactivate_output();
        gsBSpline<>& bsp = dynamic_cast<gsBSpline<>&>(*geo);
        index_t mult = 2;

        assert(mult <= bsp.degree());

        std::vector<real_t> knots;

        bsp.knots().getUniformRefinementKnots(1, knots);

        std::vector<real_t> knots2;
        knots2.reserve(mult*knots.size());

        for (size_t i = 0; i < knots.size(); ++i)
            for (int j = 0; j < mult; ++j)
                knots2.push_back(knots[i]);

        testBoehm_helper(bsp, knots2);
    }

    TEST(testCoarsening)
    {
        gsKnotVector<> kv(0.0,1.0, 7, 3,1);
        gsKnotVector<> kv_fine = kv;
        kv_fine.uniformRefine();

        gsKnotVector<> kv_coarse = kv_fine;
        std::vector<real_t> removedNodes = kv_coarse.coarsen();

        CHECK (compareKV(kv, kv_coarse));

        for (size_t i = 0; i < removedNodes.size(); ++i)
            kv_coarse.insert(removedNodes[i]);
        CHECK (compareKV(kv_coarse, kv_fine));
    }

    TEST(testTensorCoarsening)
    {
        gsKnotVector<> kv(0.0,1.0, 10,3);
        gsTensorBSplineBasis<3, real_t> tbsb(kv, kv, kv);

        gsTensorBSplineBasis<3, real_t> b_ref = tbsb;
        b_ref.uniformRefine();

        gsTensorBSplineBasis<3, real_t> b_coarse = b_ref;
        b_coarse.uniformCoarsen();

        for (unsigned i = 0; i < 3; ++i)
        {
            CHECK (compareKV(b_coarse.component(i).knots(), tbsb.component(i).knots()));
        }
    }

}
