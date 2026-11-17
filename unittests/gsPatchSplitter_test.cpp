/** @file gsPatchSplitter_test.cpp

    @brief Unit tests for gsPatchSplitter<T> (Mode A, Mode B, Polychord splitting,
           Cross-Field O-grid, and ArgyrisSplines compatibility).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Hasanova, Antigravity
**/

#include "gismo_unittest.h"
#include <gsModeling/gsPatchSplitter.h>
#include <gsModeling/gsPlanarDomainGenerator.h>

SUITE(gsPatchSplitter_test)
{
    static gsMultiPatch<real_t> makeRectMultiPatch(real_t Lx, real_t Ly)
    {
        gsKnotVector<real_t> kv(0.0, 1.0, 0, 4, 1, 3);
        gsMatrix<real_t> coefs(16, 2);
        for (int j = 0; j < 4; ++j)
            for (int i = 0; i < 4; ++i) {
                coefs(i + 4 * j, 0) = Lx * (static_cast<real_t>(i) / 3);
                coefs(i + 4 * j, 1) = Ly * (static_cast<real_t>(j) / 3);
            }
        gsMultiPatch<real_t> mp;
        mp.addPatch(gsGeometry<real_t>::uPtr(new gsTensorBSpline<2, real_t>(kv, kv, coefs)));
        mp.computeTopology();
        return mp;
    }

    static void checkStrictConformingAndPositiveDetJ(const gsMultiPatch<real_t>& mp)
    {
        CHECK(mp.nBoxes() >= 1);
        gsMultiPatch<real_t> copy(mp);
        CHECK(copy.computeTopology(1e-4, false, false));
        CHECK_EQUAL(4 * static_cast<size_t>(copy.nPatches()),
                    2 * copy.nInterfaces() + copy.nBoundary());

        for (size_t p = 0; p < copy.nPatches(); ++p) {
            gsMatrix<real_t> supp = copy.patch(p).support();
            for (int j = 0; j < 20; ++j) {
                for (int i = 0; i < 20; ++i) {
                    gsMatrix<real_t> u(2, 1), J;
                    u(0, 0) = supp(0, 0) + (supp(0, 1) - supp(0, 0)) * (i / 19.0);
                    u(1, 0) = supp(1, 0) + (supp(1, 1) - supp(1, 0)) * (j / 19.0);
                    copy.patch(p).jacobian_into(u, J);
                    real_t detJ = J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0);
                    CHECK(detJ > 1e-9);
                }
            }
        }
    }

    TEST(ModeA_StretchedRectanglePolychord)
    {
        gsMultiPatch<real_t> rect = makeRectMultiPatch(5.0, 1.0);
        gsPatchSplitter<real_t> splitter(rect);
        gsMultiPatch<real_t> res = splitter.compute();
        checkStrictConformingAndPositiveDetJ(res);
        auto qr = splitter.analyzeQuality(res);
        CHECK(qr.nPatches >= 4 && qr.nPatches <= 6);
        CHECK(qr.maxAspectRatio <= 1.3);
        CHECK(!qr.hasHangingNodes);
    }

    TEST(ModeB_UnitDiskButterflyOGrid)
    {
        std::unique_ptr<gsPlanarDomain<real_t>> disk(
            gsPlanarDomainGenerator<real_t>::preset("circle", 1.0, 3));
        gsPatchSplitter<real_t> splitter(*disk);
        splitter.options().setInt("Strategy", 1);
        gsMultiPatch<real_t> res = splitter.compute();
        checkStrictConformingAndPositiveDetJ(res);
        auto qr = splitter.analyzeQuality(res);
        CHECK_EQUAL(5, qr.nPatches);
        CHECK(qr.maxAngleDeg < 145.0);
    }

    TEST(ModeB_PlateWithHole)
    {
        std::unique_ptr<gsPlanarDomain<real_t>> dom(
            gsPlanarDomainGenerator<real_t>::preset("plate_with_hole", 1.0, 3));
        gsPatchSplitter<real_t> splitter(*dom);
        gsMultiPatch<real_t> res = splitter.compute();
        checkStrictConformingAndPositiveDetJ(res);
        auto qr = splitter.analyzeQuality(res);
        CHECK(qr.nPatches >= 4);
        CHECK(!qr.hasHangingNodes);
    }
}
