/** @file gsAdaptiveMeshingCompare_test.cpp

    @brief Tests for gsOverlapCompare::check(): a level-0 box has no
    parent, so check() must return false without ever
    reaching getParent() -- which throws GISMO_ENSURE(level()>0, ...) on a
    level-0 box (gsHSplines/gsHBox.hpp:341-346).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Verhelst
 **/

#include "gismo_unittest.h"
#include <gsHSplines/gsHBox.h>
#include <gsHSplines/gsHBoxContainer.h>
#include <gsHSplines/gsTHBSplineBasis.h>
#include <gsAssembler/gsAdaptiveMeshingCompare.h>

namespace {

// A THB basis with one level of local refinement in the lower-left corner,
// same construction pattern as unittests/gsAdaptiveParametrizatin_test.cpp
// F2 (:602-606): kv(0,1,3,3) is degree 2 with 3 interior knots (4 elements
// per direction at level 0); refineElements({{1,0,0,2,2}}) refines the
// level-1 box [0,2]^2 == parameter [0,0.25]^2, creating one level of local
// refinement there.
gsTHBSplineBasis<2,real_t> makeRefinedBasis()
{
    gsKnotVector<real_t> kv(0, 1, 3, 3);
    gsTensorBSplineBasis<2,real_t> tb(kv, kv);
    gsTHBSplineBasis<2,real_t> thb(tb);
    thb.refineElements({{1, 0, 0, 2, 2}});
    return thb;
}

} // anonymous namespace

SUITE(gsAdaptiveMeshingCompare_test)
{

    TEST(OverlapCompare_Level0BoxReturnsFalse)
    {
        gsTHBSplineBasis<2,real_t> thb = makeRefinedBasis();

        // A genuine level-0 leaf away from the refined corner.
        gsVector<index_t,2> low, upp;
        low << 2, 2; upp << 3, 3;
        gsHBox<2,real_t> level0Box(low, upp, 0, &thb);

        gsHBoxContainer<2,real_t> markedRef(level0Box);
        gsOverlapCompare<2,real_t> cmp(markedRef, /*m*/2);

        // Without the guard (box.level()==0 return false;), check() would
        // reach getParent() on a level-0 box and throw std::runtime_error
        // (GISMO_ENSURE(this->level()>0, ...), gsHBox.hpp:341-346).
        // UnitTest++ fails a test on an escaping exception, so an
        // exception here is itself a falsification signal, but the guard's
        // documented behaviour is a returned false, which is what this
        // CHECK_EQUAL pins.
        CHECK_EQUAL(false, cmp.check(level0Box));
    }

    TEST(OverlapCompare_CtorAcceptsLevel0MarkedBox)
    {
        gsTHBSplineBasis<2,real_t> thb = makeRefinedBasis();

        gsVector<index_t,2> low, upp;
        low << 2, 2; upp << 3, 3;
        gsHBox<2,real_t> level0Box(low, upp, 0, &thb);
        gsHBoxContainer<2,real_t> markedRef(level0Box);

        // The CONSTRUCTOR only calls gsHBoxContainer::getChildren(), which
        // delegates to gsHBox::getChildren() (gsHBox.hpp:372-377) -- that
        // has no level restriction, unlike getParent() (the sibling call
        // the check() guard protects). Constructing from a level-0 marked
        // box must not throw.
        gsOverlapCompare<2,real_t> cmp(markedRef, 2);
        CHECK(true);   // reaching here means construction did not throw
    }

    TEST(OverlapCompare_Level1BoxReachesParentPath)
    {
        gsTHBSplineBasis<2,real_t> thb = makeRefinedBasis();

        gsVector<index_t,2> low0, upp0;
        low0 << 2, 2; upp0 << 3, 3;
        gsHBox<2,real_t> level0Box(low0, upp0, 0, &thb);
        gsHBoxContainer<2,real_t> markedRef(level0Box);
        gsOverlapCompare<2,real_t> cmp(markedRef, 2);

        // A level-1 leaf inside the refined corner.
        gsVector<index_t,2> low1, upp1;
        low1 << 0, 0; upp1 << 1, 1;
        gsHBox<2,real_t> level1Box(low1, upp1, 1, &thb);

        // check() on a level-1 box bypasses the guard and reaches
        // getParent()/getCextension()/levelInCenter(). The returned value
        // depends on those in a way not obvious from this fixture, so it is
        // OBSERVED here, not asserted -- all the discriminating power of
        // guard 3 lives in the level-0 test above and its falsification.
        // The only thing this test pins is that the getParent() path does
        // not throw for a legitimate level-1 box (an uncaught exception
        // here would still fail the test).
        bool result = cmp.check(level1Box);
        gsTestInfo << "OverlapCompare_Level1BoxReachesParentPath: check() = "
                   << (result ? "true" : "false") << "\n";
    }

}
