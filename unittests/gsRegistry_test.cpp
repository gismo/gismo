/** @file gsRegistry_test.cpp

    @brief Tests for the generic gsRegistry name-to-factory registry

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "gismo_unittest.h"
#include <gsCore/gsRegistry.h>
#include <gsCore/gsRegistry.hpp> // local instantiation for the test base

namespace
{
struct TestBase
{
    virtual ~TestBase() { }
    virtual int id() const = 0;
    typedef gismo::memory::unique_ptr<TestBase> uPtr;
};
struct TestA : public TestBase { int id() const { return 1; } };
struct TestB : public TestBase { int id() const { return 2; } };
}

SUITE(gsRegistry_test)
{

TEST(register_make_list)
{
    gsRegistry<TestBase> & reg = gsRegistry<TestBase>::get();

    CHECK( !reg.has("A") );
    reg.add("A", gsRegistryFactory<TestBase,TestA>);
    reg.add("B", gsRegistryFactory<TestBase,TestB>);
    CHECK( reg.has("A") );
    CHECK( reg.has("B") );

    CHECK_EQUAL( 1, reg.make("A")->id() );
    CHECK_EQUAL( 2, reg.make("B")->id() );
    CHECK_EQUAL( 2, (int)reg.list().size() );

    // same instance on re-access
    CHECK( &reg == &gsRegistry<TestBase>::get() );
}

TEST(idempotent_registration)
{
    gsRegistry<TestBase> & reg = gsRegistry<TestBase>::get();
    // re-registering "A" with another factory is a no-op
    reg.add("A", gsRegistryFactory<TestBase,TestB>);
    CHECK_EQUAL( 1, reg.make("A")->id() );
}

TEST(missing_name_throws)
{
    CHECK_THROW( gsRegistry<TestBase>::get().make("does-not-exist"),
                 std::exception );
}

}
