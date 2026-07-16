/** @file gsParaview_test.cpp

    @brief Tests for gsParaview class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include "gismo_unittest.h"
#include <gsIO/gsParaview.h>

SUITE(gsParaview_test)
{

TEST(DefaultOptions)
{
    gsParaview<real_t> pv;
    CHECK_EQUAL(1000, pv.options().getInt("numPoints"));
    CHECK_EQUAL(12,   pv.options().getInt("precision"));
    CHECK(!pv.options().getSwitch("plotElements"));
    CHECK(!pv.options().getSwitch("plotControlNet"));
    CHECK(!pv.options().getSwitch("show"));
}

TEST(WriteGeometry_smoke)
{
    const std::string tmp = gsFileManager::getTempPath();
    if (tmp.empty()) return;

    const std::string fn = tmp + "gsParaview_geo_test";

    gsMultiPatch<> mp;
    mp.addPatch(gsNurbsCreator<>::BSplineSquare());

    gsParaview<real_t> pv;
    pv.write(mp, fn);

    CHECK(gsFileManager::fileExists(fn + ".pvd"));
}

TEST(WriteMultiPatch_options_smoke)
{
    const std::string tmp = gsFileManager::getTempPath();
    if (tmp.empty()) return;

    const std::string fn = tmp + "gsParaview_mp_test";

    gsMultiPatch<> mp;
    mp.addPatch(gsNurbsCreator<>::BSplineSquare());
    mp.computeTopology();

    gsParaview<real_t> pv;
    pv.options().setSwitch("plotElements",   true);
    pv.options().setSwitch("plotControlNet", true);
    pv.write(mp, fn);

    CHECK(gsFileManager::fileExists(fn + ".pvd"));
    // mesh and control net files should also be produced
    CHECK(gsFileManager::fileExists(fn + "_0_mesh.vtp") ||
          gsFileManager::fileExists(fn + "_patch0_mesh.vtp") ||
          gsFileManager::fileExists(fn + "_0.vts"));
}

TEST(PrecisionOption_honored)
{
    const std::string tmp = gsFileManager::getTempPath();
    if (tmp.empty()) return;

    gsMultiPatch<> mp;
    mp.addPatch(gsNurbsCreator<>::BSplineSquare());

    const std::string fn3  = tmp + "gsParaview_prec3";
    const std::string fn12 = tmp + "gsParaview_prec12";

    {
        gsParaview<real_t> pv;
        pv.options().setInt("precision", 3);
        pv.write(mp, fn3);
    }
    {
        gsParaview<real_t> pv;
        pv.options().setInt("precision", 12);
        pv.write(mp, fn12);
    }

    // Higher precision produces a strictly larger file
    const std::string vts3  = fn3  + "0.vts";
    const std::string vts12 = fn12 + "0.vts";

    if (gsFileManager::fileExists(vts3) && gsFileManager::fileExists(vts12))
    {
        std::ifstream f3(vts3);
        std::ifstream f12(vts12);
        std::string s3((std::istreambuf_iterator<char>(f3)), std::istreambuf_iterator<char>());
        std::string s12((std::istreambuf_iterator<char>(f12)), std::istreambuf_iterator<char>());
        CHECK(s12.size() > s3.size());
    }
}

TEST(WriteField_smoke)
{
    const std::string tmp = gsFileManager::getTempPath();
    if (tmp.empty()) return;

    const std::string fn = tmp + "gsParaview_field_test";

    gsMultiPatch<> mp;
    mp.addPatch(gsNurbsCreator<>::BSplineSquare());

    // Identity field: function = geometry
    gsField<> field(mp, mp);

    gsParaview<real_t> pv;
    pv.write(field, fn);

    CHECK(gsFileManager::fileExists(fn + ".pvd"));
}

TEST(TimeSteppingWorkflow_smoke)
{
    const std::string tmp = gsFileManager::getTempPath();
    if (tmp.empty()) return;

    const std::string fn = tmp + "gsParaview_ts_test";

    gsMultiPatch<> mp;
    mp.addPatch(gsNurbsCreator<>::BSplineSquare());

    // Identity field: function = geometry
    gsField<> field(mp, mp);

    gsParaviewCollection collection(fn);
    collection.newTimeStep(&mp, 0.0);
    collection.addField(field, "identity");
    collection.saveTimeStep();
    collection.save();

    CHECK(gsFileManager::fileExists(fn + ".pvd"));
}

} // SUITE
