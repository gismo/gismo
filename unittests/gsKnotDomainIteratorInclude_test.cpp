/** @file gsKnotDomainIteratorInclude_test.cpp

    @brief Checks that gsKnotDomainIterator.h compiles when included first

    The include of <gsDomain/gsKnotDomainIterator.h> must stay the first
    include of this file: the suite verifies that the header compiles when it
    is included before any other gismo header.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Verhelst
**/

#include <gsDomain/gsKnotDomainIterator.h>

#include "gismo_unittest.h"

SUITE(gsKnotDomainIteratorInclude_test)
{
    TEST(iterateKnotVectorElements)
    {
        // Degree 2, knots 0 0 0 .25 .5 .75 1 1 1
        gsKnotVector<real_t> kv(0, 1, 3, 3);

        gsKnotVector<real_t>::domainIter it    = kv.beginAll();
        gsKnotVector<real_t>::domainIter itEnd = kv.endAll();
        const std::vector<real_t> br = kv.breaks();

        size_t count = 0;
        for (; it != itEnd; ++it)
        {
            CHECK( count + 1 < br.size() );
            if (count + 1 >= br.size()) break;
            CHECK_CLOSE( br[count],   it.lowerCorner()(0), 1e-12 );
            CHECK_CLOSE( br[count+1], it.upperCorner()(0), 1e-12 );
            ++count;
        }

        CHECK_EQUAL( kv.numElements(), count );
        CHECK_EQUAL( size_t(4), count );
    }
}
