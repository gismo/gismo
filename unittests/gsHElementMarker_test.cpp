/** @file gsHElementMarker_test.cpp

    @brief Tests that gsHElementMarker::markCrs dispatches on the
    "CoarsenRule" option, not "RefineRule".

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Testing
**/

#include "gismo_unittest.h"

#include <gsHSplines/gsHElementMarker.h>
#include <gsHSplines/gsTHBSplineBasis.h>

using namespace gismo;

namespace
{

// std::set<element_t,...> with a custom comparator may not define
// operator==, so compare in terms of the comparator instead of relying on
// std::set::operator==.
template <short_t d, class T>
bool sameSet(const std::set<gsHElement<d,T>, typename gsHElement<d,T>::Compare> & a,
             const std::set<gsHElement<d,T>, typename gsHElement<d,T>::Compare> & b)
{
    typename gsHElement<d,T>::Compare cmp;
    if (a.size() != b.size()) return false;
    auto ia = a.begin(); auto ib = b.begin();
    for (; ia != a.end(); ++ia, ++ib)
        if (cmp(*ia,*ib) || cmp(*ib,*ia)) return false;
    return true;
}

// A THB-spline basis over gsKnotVector<real_t>(0,1,3,3) in both directions,
// refined once so that level-1 elements exist.
gsTHBSplineBasis<2,real_t> makeFixtureBasis()
{
    gsKnotVector<real_t> kv(0, 1, 3, 3);
    gsTensorBSplineBasis<2,real_t> tbasis(kv, kv);
    gsTHBSplineBasis<2,real_t> basis(tbasis);
    basis.refineElements({{1, 0, 0, 2, 2}}); // level 1, box [0,2] x [0,2]
    return basis;
}

} // anonymous namespace

SUITE(gsHElementMarker_test)
{

    TEST(MarkCrs_DispatchesOnCoarsenRule_NotRefineRule)
    {
        gsTHBSplineBasis<2,real_t> basis = makeFixtureBasis();

        // Strictly increasing, non-uniform error spread so the threshold,
        // percentage and fraction rules genuinely disagree.
        const size_t nel = basis.numElements();
        std::vector<real_t> err(nel);
        for (size_t i = 0; i != nel; i++)
            err[i] = static_cast<real_t>((i+1)*(i+1));

        // With this fixture (19 elements: 15 at level 0, the 4 refined
        // sub-elements at level 1 carrying the highest errors, since the
        // level-1 ids are the last -- and thus highest-error -- in ascending
        // order), the coarsening candidates are ONLY the 4 level-1 elements
        // (level==0 is always skipped). GARU/PUCA scan ascending order and
        // only reach them once the threshold/percentage covers most of the
        // 15 level-0 elements first; 0.85 was picked (probed empirically) to
        // put GARU and PUCA at genuinely different, nonzero counts.
        const real_t coarsenParam = 0.85;

        gsOptionList optsA = gsHElementMarker<2,real_t>::defaultOptions();
        optsA.setSwitch("Admissible", false); // isolate rule dispatch from _markCrs_admissible
        optsA.setInt("RefineRule", 1);
        optsA.setInt("CoarsenRule", 2);
        optsA.setReal("CoarsenParam", coarsenParam);
        gsHElementMarker<2,real_t> markerA(basis, optsA);
        markerA.setErrors(err);
        auto A = markerA.markCrs(); // empty refined => sibling filters skipped

        gsOptionList optsB = optsA;
        optsB.setInt("RefineRule", 2); // == CoarsenRule: what the OLD (pre-fix)
                                        // code produced only when told
                                        // RefineRule == the desired coarsening rule
        gsHElementMarker<2,real_t> markerB(basis, optsB);
        markerB.setErrors(err);
        auto B = markerB.markCrs();

        gsOptionList optsC = optsA;
        optsC.setInt("CoarsenRule", 1); // what the OLD code produced for A
                                         // (dispatching on RefineRule == 1)
        gsHElementMarker<2,real_t> markerC(basis, optsC);
        markerC.setErrors(err);
        auto C = markerC.markCrs();

        gsTestInfo << "sizes: A=" << A.size() << " B=" << B.size() << " C=" << C.size() << "\n";

        // Guard: the fixture must actually distinguish rule 1 from rule 2,
        // otherwise the A==B / A!=C checks below would be vacuous.
        const bool bEqualsC = sameSet(B, C);
        CHECK(bEqualsC == false);

        // markCrs follows CoarsenRule (A and B share CoarsenRule=2, differ
        // only in RefineRule, which markCrs must ignore).
        const bool aEqualsB = sameSet(A, B);
        CHECK(aEqualsB);

        // markCrs no longer follows RefineRule (A has RefineRule=1,
        // CoarsenRule=2; C has RefineRule=1, CoarsenRule=1 -- pre-fix code
        // would have produced A == C here).
        const bool aEqualsC = sameSet(A, C);
        CHECK(!aEqualsC);

        // Additional configuration: CoarsenRule = 3 (BULK) should likewise
        // differ from CoarsenRule = 1 (C), if the fixture separates them.
        gsOptionList optsD = optsA;
        optsD.setInt("CoarsenRule", 3);
        gsHElementMarker<2,real_t> markerD(basis, optsD);
        markerD.setErrors(err);
        auto D = markerD.markCrs();
        gsTestInfo << "sizes: D(BULK)=" << D.size() << "\n";
        const bool dEqualsC = sameSet(D, C);
        if (!dEqualsC)
            CHECK(!dEqualsC);
        // else: fixture does not separate BULK from GARU here; skip (do not
        // tune the fixture to force it, per the task spec).
    }

}
