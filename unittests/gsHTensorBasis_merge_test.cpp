/** @file gsHTensorBasis_merge_test.cpp

    @brief Tests for gsHTensorBasis::merge (both instance and static factory).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Testing
**/

#include "gismo_unittest.h"

using namespace gismo;

SUITE(gsHTensorBasis_merge_test)
{
    // Create a simple 2D THB-spline basis for testing
    gsTHBSplineBasis<2>* makeTestBasis()
    {
        // Create a level-0 tensor-product basis on [0,1]^2
        // Degree 2, 4 elements per direction (5 unique knots per direction)
        std::vector<gsTensorBSplineBasis<2, double>> bases;
        gsKnotVector<double> knots = gsKnotVector<double>(0, 1, 5, 2);
        bases.push_back(gsTensorBSplineBasis<2, double>(knots, knots));

        return new gsTHBSplineBasis<2>(bases[0]);
    }

    TEST(MergeWithSelf)
    {
        // Merging a basis with itself should preserve properties
        gsTHBSplineBasis<2>* basis = makeTestBasis();

        index_t orig_size = basis->size();
        unsigned orig_maxLevel = basis->tree().getMaxInsLevel();

        // Merge with itself
        basis->merge(*basis);

        CHECK_EQUAL(basis->size(), orig_size);
        CHECK_EQUAL(basis->tree().getMaxInsLevel(), orig_maxLevel);

        delete basis;
    }

    TEST(MergeIsAtLeastAsRefined)
    {
        // After merging two bases, the result should be at least as refined as either input
        gsTHBSplineBasis<2>* basis1 = makeTestBasis();
        gsTHBSplineBasis<2>* basis2 = makeTestBasis();

        index_t size1 = basis1->size();
        index_t size2 = basis2->size();
        unsigned maxLevel1 = basis1->tree().getMaxInsLevel();

        // Refine basis2 in a region
        std::array<index_t, 2> lo = {0, 0};
        std::array<index_t, 2> up = {2, 2};
        basis2->refineElements({{1, lo[0], lo[1], up[0], up[1]}}); // level 1, box [0,2] x [0,2]

        unsigned maxLevel2 = basis2->tree().getMaxInsLevel();

        // Now merge
        basis1->merge(*basis2);

        CHECK(basis1->size() >= size1);
        CHECK(basis1->size() >= size2);
        CHECK(basis1->tree().getMaxInsLevel() >= std::max(maxLevel1, maxLevel2));

        delete basis1;
        delete basis2;
    }

    TEST(MergeIsCommutative)
    {
        // Merging A with B should give same size/maxLevel as merging B with A
        gsTHBSplineBasis<2>* basisA = makeTestBasis();
        gsTHBSplineBasis<2>* basisB = makeTestBasis();

        // Refine A
        std::array<index_t, 2> loA = {0, 0};
        std::array<index_t, 2> upA = {2, 2};
        basisA->refineElements({{1, loA[0], loA[1], upA[0], upA[1]}});

        // Refine B differently
        std::array<index_t, 2> loB = {2, 2};
        std::array<index_t, 2> upB = {4, 4};
        basisB->refineElements({{1, loB[0], loB[1], upB[0], upB[1]}});

        // Store original states
        gsTHBSplineBasis<2>* origA = basisA->clone().release();
        gsTHBSplineBasis<2>* origB = basisB->clone().release();

        // Merge A + B
        basisA->merge(*basisB);
        index_t sizeAB = basisA->size();
        unsigned maxLevelAB = basisA->tree().getMaxInsLevel();

        // Merge B + A
        origB->merge(*origA);
        index_t sizeBA = origB->size();
        unsigned maxLevelBA = origB->tree().getMaxInsLevel();

        CHECK_EQUAL(sizeAB, sizeBA);
        CHECK_EQUAL(maxLevelAB, maxLevelBA);

        delete basisA;
        delete basisB;
        delete origA;
        delete origB;
    }

    TEST(StaticMergeMatchesMutatingMerge)
    {
        // The static factory merge should produce the same result as the instance method
        gsTHBSplineBasis<2>* basis1 = makeTestBasis();
        gsTHBSplineBasis<2>* basis2 = makeTestBasis();

        // Refine both
        std::array<index_t, 2> lo = {0, 0};
        std::array<index_t, 2> up = {2, 2};
        basis1->refineElements({{1, lo[0], lo[1], up[0], up[1]}});
        basis2->refineElements({{1, lo[0], lo[1], up[0], up[1]}});

        // Store copies for static call
        gsTHBSplineBasis<2>* staticBasis1 = basis1->clone().release();
        gsTHBSplineBasis<2>* staticBasis2 = basis2->clone().release();

        // Instance method
        basis1->merge(*basis2);
        index_t mutating_size = basis1->size();
        unsigned mutating_maxLevel = basis1->tree().getMaxInsLevel();

        // Static method
        auto static_result = gsTHBSplineBasis<2>::merge(*staticBasis1, *staticBasis2);
        index_t static_size = static_result->size();
        unsigned static_maxLevel = static_result->tree().getMaxInsLevel();

        CHECK_EQUAL(mutating_size, static_size);
        CHECK_EQUAL(mutating_maxLevel, static_maxLevel);

        delete basis1;
        delete basis2;
        delete staticBasis1;
        delete staticBasis2;
    }

    TEST(MergedBasisPartitionOfUnity)
    {
        // After merge, the basis should still form a partition of unity
        gsTHBSplineBasis<2>* basis = makeTestBasis();

        // Refine to create non-trivial hierarchy
        std::array<index_t, 2> lo = {0, 0};
        std::array<index_t, 2> up = {2, 2};
        basis->refineElements({{1, lo[0], lo[1], up[0], up[1]}});

        gsTHBSplineBasis<2>* basis2 = basis->clone().release();
        std::array<index_t, 2> lo2 = {2, 2};
        std::array<index_t, 2> up2 = {4, 4};
        basis2->refineElements({{1, lo2[0], lo2[1], up2[0], up2[1]}});

        // Merge
        basis->merge(*basis2);

        // Test partition of unity: should pass at level 1e-10
        CHECK(basis->testPartitionOfUnity(100, 1e-10));

        delete basis;
        delete basis2;
    }
}
