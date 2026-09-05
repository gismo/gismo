/** @file gsXmlRegistry_test.cpp

    @brief Tests for the XML reader/writer registry primitives

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "gismo_unittest.h"
#include <gsIO/gsXmlRegistry.h>

using gismo::internal::gsXmlRegistry;

namespace
{
void fnA() { }
void fnB() { }
void fnC() { }
}

SUITE(gsXmlRegistry_test)
{

TEST(get_registration)
{
    gsXmlRegistry & r = gsXmlRegistry::get();
    CHECK( &r == &gsXmlRegistry::get() );

    CHECK( NULL == r.findGet("TestBase", "TypeA") );
    r.addGet("TestBase", "TypeA", &fnA);
    r.addGet("TestBase", "TypeB", &fnB);
    CHECK( &fnA == r.findGet("TestBase", "TypeA") );
    CHECK( &fnB == r.findGet("TestBase", "TypeB") );
    // idempotent: first registration wins
    r.addGet("TestBase", "TypeA", &fnB);
    CHECK( &fnA == r.findGet("TestBase", "TypeA") );
    // base namespaces are separate
    CHECK( NULL == r.findGet("OtherBase", "TypeA") );
}

TEST(put_chain_priorities)
{
    gsXmlRegistry & r = gsXmlRegistry::get();
    // register out of order; chain must come back priority-sorted
    r.addPut("TestBase2", 200, &fnB);
    r.addPut("TestBase2", 100, &fnA);
    r.addPut("TestBase2", 900, &fnC);
    r.addPut("TestBase2", 100, &fnA); // idempotent

    std::vector<gsXmlRegistry::AnyFn> chain;
    r.putChain("TestBase2", chain);
    CHECK_EQUAL( 3, (int)chain.size() );
    CHECK( &fnA == chain[0] );
    CHECK( &fnB == chain[1] );
    CHECK( &fnC == chain[2] );
}

}
