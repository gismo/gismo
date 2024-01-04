/** @file gsFixedIssues_test.cpp

    @brief Unit tests for solved issues on github.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: D. Mokris

**/

#include "gismo_unittest.h"

using namespace gismo;

SUITE(gsFixedIssues_test)
{
	TEST(issue_672)
	{
		typedef index_t domain_index_type;
		gsKnotVector<real_t> uKnots(0, 1, 127, 3);
		gsKnotVector<real_t> vKnots(0, 1, 15, 3);
		gsTensorBSplineBasis<2, real_t> tensBasis(uKnots, vKnots);
		gsTHBSplineBasis<2> thbBasis(tensBasis);

		unsigned indexLevel = thbBasis.tree().getIndexLevel();
		gsVector<domain_index_type, 2> upp = thbBasis.tree().upperCorner();
		domain_index_type max = std::numeric_limits<domain_index_type>::max();

		// Check upp against under- and overflow
		CHECK(upp[0] >  0);
		CHECK(upp[0] <= max);

		CHECK(upp[1] >  0);
		CHECK(upp[1] <= max);

		// Until 677 is fixed, check that there is still place for the extra local2global.
		CHECK((upp[0] << indexLevel) <= max);
		CHECK((upp[1] << indexLevel) <= max);
	}
}
