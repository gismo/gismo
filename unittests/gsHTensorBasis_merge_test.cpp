/** @file gsHTensorBasis_merge_test.cpp

    @brief Tests for gsHTensorBasis::merge (both instance and static factory).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Testing
**/

#include "gismo_unittest.h"

#include <algorithm>
#include <array>

using namespace gismo;

namespace
{
    struct BoxSignature
    {
        index_t level;
        std::array<index_t, 2> lower;
        std::array<index_t, 2> upper;

        bool operator==(const BoxSignature &other) const
        {
            return level == other.level && lower == other.lower && upper == other.upper;
        }

        bool operator<(const BoxSignature &other) const
        {
            if (level != other.level)
                return level < other.level;
            if (lower != other.lower)
                return lower < other.lower;
            return upper < other.upper;
        }
    };

    template <class Basis>
    std::vector<BoxSignature> boxSignatures(const Basis &basis)
    {
        gsMatrix<index_t> lower, upper;
        gsVector<index_t> level;
        basis.tree().getBoxes(lower, upper, level);

        std::vector<BoxSignature> result;
        result.reserve(static_cast<size_t>(level.size()));
        for (index_t i = 0; i < level.size(); ++i)
        {
            BoxSignature sig;
            sig.level = level[i];
            for (short_t j = 0; j < 2; ++j)
            {
                sig.lower[static_cast<size_t>(j)] = lower(i, j);
                sig.upper[static_cast<size_t>(j)] = upper(i, j);
            }
            result.push_back(sig);
        }

        std::sort(result.begin(), result.end());
        return result;
    }
}

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
        // Merging a basis with itself should preserve the topology.
        gsTHBSplineBasis<2>* basis = makeTestBasis();

        const index_t orig_size = basis->size();
        const unsigned orig_maxLevel = basis->tree().getMaxInsLevel();
        const std::vector<BoxSignature> original = boxSignatures(*basis);

        // Merge with itself
        basis->merge(*basis);

        CHECK_EQUAL(basis->size(), orig_size);
        CHECK_EQUAL(basis->tree().getMaxInsLevel(), orig_maxLevel);
        CHECK(boxSignatures(*basis) == original);

        delete basis;
    }

    TEST(MergeIsAtLeastAsRefined)
    {
        // Merging should recover the finer basis when the other input is unrefined.
        gsTHBSplineBasis<2>* basis1 = makeTestBasis();
        gsTHBSplineBasis<2>* basis2 = makeTestBasis();

        // Refine basis2 in a region
        std::array<index_t, 2> lo = {0, 0};
        std::array<index_t, 2> up = {2, 2};
        basis2->refineElements({{1, lo[0], lo[1], up[0], up[1]}}); // level 1, box [0,2] x [0,2]

        const index_t size2 = basis2->size();
        const unsigned maxLevel2 = basis2->tree().getMaxInsLevel();

        // Now merge
        basis1->merge(*basis2);

        CHECK_EQUAL(basis1->size(), size2);
        CHECK_EQUAL(basis1->tree().getMaxInsLevel(), maxLevel2);
        CHECK(boxSignatures(*basis1) == boxSignatures(*basis2));

        delete basis1;
        delete basis2;
    }

    TEST(MergeIsCommutative)
    {
        // Merging A with B should give the same merged topology as B with A.
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

        gsTHBSplineBasis<2>* expected = makeTestBasis();
        expected->refineElements({{1, loA[0], loA[1], upA[0], upA[1]}});
        expected->refineElements({{1, loB[0], loB[1], upB[0], upB[1]}});

        // Merge A + B
        basisA->merge(*basisB);
        const index_t sizeAB = basisA->size();
        const unsigned maxLevelAB = basisA->tree().getMaxInsLevel();
        const std::vector<BoxSignature> sigAB = boxSignatures(*basisA);

        // Merge B + A
        origB->merge(*origA);
        const index_t sizeBA = origB->size();
        const unsigned maxLevelBA = origB->tree().getMaxInsLevel();
        const std::vector<BoxSignature> sigBA = boxSignatures(*origB);

        CHECK_EQUAL(sizeAB, sizeBA);
        CHECK_EQUAL(maxLevelAB, maxLevelBA);
        CHECK(sigAB == sigBA);
        CHECK(sigAB == boxSignatures(*expected));

        delete basisA;
        delete basisB;
        delete origA;
        delete origB;
        delete expected;
    }

    TEST(StaticMergeMatchesMutatingMerge)
    {
        // The static factory merge should produce the same topology as the instance method.
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
        const index_t mutating_size = basis1->size();
        const unsigned mutating_maxLevel = basis1->tree().getMaxInsLevel();
        const std::vector<BoxSignature> mutating_sig = boxSignatures(*basis1);

        // Static method
        auto static_result = gsTHBSplineBasis<2>::merge(*staticBasis1, *staticBasis2);
        const index_t static_size = static_result->size();
        const unsigned static_maxLevel = static_result->tree().getMaxInsLevel();
        const std::vector<BoxSignature> static_sig = boxSignatures(*static_result);

        CHECK_EQUAL(mutating_size, static_size);
        CHECK_EQUAL(mutating_maxLevel, static_maxLevel);
        CHECK(mutating_sig == static_sig);

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

        // Test partition of unity: should pass at level 1e-10.
        CHECK(basis->testPartitionOfUnity(100, 1e-10));

        delete basis;
        delete basis2;
    }
}
