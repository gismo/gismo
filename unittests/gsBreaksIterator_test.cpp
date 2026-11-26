/** @file gsBreaksIterator_test.cpp

    @brief Tests for gsBreaksIterator: boundary-element detection and clamping of strides past the end

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Verhelst
*/

#include "gismo_unittest.h"
#include <gsDomain/gsBreaksIterator.h>

using namespace gismo;

SUITE(gsBreaksIterator_test)
{
    TEST(IsBoundaryElement_MultiElement)
    {
        std::vector<real_t> knots; // 4 elements
        knots.push_back(0.0); knots.push_back(0.25);
        knots.push_back(0.5); knots.push_back(0.75); knots.push_back(1.0);

        gsDomainIteratorWrapper<real_t> endIt = gsBreaksIterator<real_t>::make(knots, false);
        std::vector<bool> got;
        for (gsDomainIteratorWrapper<real_t> it = gsBreaksIterator<real_t>::make(knots, true);
             it != endIt; ++it)
            got.push_back(it.isBoundaryElement());

        CHECK ( got.size() == 4 );
        CHECK (  got[0] );   // first element
        CHECK ( !got[1] );   // interior
        CHECK ( !got[2] );   // interior
        CHECK (  got[3] );   // last element
    }

    TEST(IsBoundaryElement_SingleElement)
    {
        std::vector<real_t> knots; // 1 element: first == last
        knots.push_back(0.0); knots.push_back(1.0);

        gsDomainIteratorWrapper<real_t> it = gsBreaksIterator<real_t>::make(knots, true);
        CHECK ( it.isBoundaryElement() );
    }

    // A stride that overshoots must stop at the last break (the end
    // position) instead of running past it. There are no entries beyond the
    // last break to read, so the stop position is observed by stepping back
    // once: from the end position the iterator then sits on the last element,
    // [0.75, 1], with lowerCorner 0.75 (the fourth break, copied unchanged
    // from the input, hence an exact comparison). Had the stride run one
    // break past the end, stepping back would land on the end position, whose
    // lowerCorner is 1 and which is not a boundary element. The element id is
    // not checked: it counts the steps taken and is not clamped.
    TEST(StridePastEndStopsAtEnd)
    {
        std::vector<real_t> knots; // 4 elements
        knots.push_back(0.0); knots.push_back(0.25);
        knots.push_back(0.5); knots.push_back(0.75); knots.push_back(1.0);

        gsDomainIteratorWrapper<real_t> it = gsBreaksIterator<real_t>::make(knots, true);
        it += 3; // last element
        it += 2; // overshoots the end position by one
        --it;    // back onto the last element if the stride was clamped

        CHECK_EQUAL ( knots[3], it.lowerCorner()(0) );
        CHECK ( it.isBoundaryElement() );
    }
}
