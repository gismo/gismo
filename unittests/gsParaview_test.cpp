/** @file gsParaview_test.cpp

    @brief Tests for gsParaview class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include "gismo_unittest.h"
#include <gsIO/gsParaviewCollection.h>

#include <cstdio> // std::remove
#include <sstream> // std::istringstream

SUITE(gsParaview_test)
{

TEST(DefaultOptions)
{
    gsParaview<real_t> pv;
    CHECK_EQUAL(1000, pv.options().getInt("numPoints"));
    CHECK_EQUAL(5,    pv.options().getInt("precision"));
    CHECK(!pv.options().getSwitch("elements"));
    CHECK(!pv.options().getSwitch("controlNet"));
    CHECK(!pv.options().getSwitch("singleFile"));
    CHECK(!pv.options().getSwitch("base64"));
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
    pv.options().setSwitch("elements",   true);
    pv.options().setSwitch("controlNet", true);
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

TEST(WriteSingleFileVTU_smoke)
{
    const std::string tmp = gsFileManager::getTempPath();
    if (tmp.empty()) return;

    gsMultiPatch<> mp;
    mp.addPatch(gsNurbsCreator<>::BSplineSquare());
    mp.addPatch(gsNurbsCreator<>::BSplineSquare());

    {
        const std::string fn = tmp + "gsParaview_mp_singlefile_test";
        gsParaview<real_t> pv;
        pv.options().setSwitch("singleFile", true);
        pv.write(mp, fn);

        CHECK(gsFileManager::fileExists(fn + ".vtu"));
        std::ifstream f(fn + ".vtu");
        std::string s((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
        CHECK(s.find("<UnstructuredGrid>") != std::string::npos);
        CHECK(s.find("PatchID") != std::string::npos);
    }

    {
        const std::string fn = tmp + "gsParaview_field_singlefile_test";
        gsField<> field(mp, mp);
        gsParaview<real_t> pv;
        pv.options().setSwitch("singleFile", true);
        pv.write(field, fn);

        CHECK(gsFileManager::fileExists(fn + ".vtu"));
        std::ifstream f(fn + ".vtu");
        std::string s((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
        CHECK(s.find("<UnstructuredGrid>") != std::string::npos);
        CHECK(s.find("SolutionField") != std::string::npos);
    }
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

    gsParaviewCollection<real_t> collection(fn);
    collection.newTimeStep(mp, 0.0);
    collection.addField(field, "identity");
    collection.saveTimeStep();
    collection.save();

    CHECK(gsFileManager::fileExists(fn + ".pvd"));
}

TEST(SingleFileExtras_mergedAcrossPatches)
{
    // With "singleFile", the element mesh and the control net are merged over
    // all patches into one file each, so the output count does not grow with
    // nPatches(). Two patches must still give exactly one _mesh and one _cnet.
    const std::string tmp = gsFileManager::getTempPath();
    if (tmp.empty()) return;

    const std::string fn = tmp + "gsParaview_sf_merged_test";

    // getTempPath() falls back to the working directory when TMPDIR is unset,
    // so outputs persist between runs: the negative checks below would then
    // report on a leftover file rather than on this run.
    std::remove((fn + ".pvd").c_str());
    std::remove((fn + ".vtu").c_str());
    std::remove((fn + "_mesh.vtp").c_str());
    std::remove((fn + "_cnet.vtp").c_str());
    std::remove((fn + "_0_mesh.vtp").c_str());
    std::remove((fn + "_1_mesh.vtp").c_str());

    gsMultiPatch<> mp;
    mp.addPatch(gsNurbsCreator<>::BSplineSquare());
    mp.addPatch(gsNurbsCreator<>::BSplineSquare(1, 1, 0));
    mp.computeTopology();

    gsParaview<real_t> pv;
    pv.options().setSwitch("singleFile", true);
    pv.options().setSwitch("elements",   true);
    pv.options().setSwitch("controlNet", true);
    pv.write(mp, fn);

    CHECK(gsFileManager::fileExists(fn + ".pvd"));
    CHECK(gsFileManager::fileExists(fn + ".vtu"));
    CHECK(gsFileManager::fileExists(fn + "_mesh.vtp"));
    CHECK(gsFileManager::fileExists(fn + "_cnet.vtp"));
    // the per-patch files of the non-merged path must NOT appear
    CHECK(!gsFileManager::fileExists(fn + "_0_mesh.vtp"));
    CHECK(!gsFileManager::fileExists(fn + "_1_mesh.vtp"));
}

TEST(SingleFileExtras_noPvd)
{
    // "writePvd"=false together with the mesh/control-net switches: the two
    // extra files must still be written, and no collection file created.
    const std::string tmp = gsFileManager::getTempPath();
    if (tmp.empty()) return;

    const std::string fn = tmp + "gsParaview_sf_nopvd_test";

    // A stale .pvd from an earlier run would make the negative check below
    // pass or fail for the wrong reason.
    std::remove((fn + ".pvd").c_str());
    std::remove((fn + ".vtu").c_str());
    std::remove((fn + "_mesh.vtp").c_str());
    std::remove((fn + "_cnet.vtp").c_str());

    gsMultiPatch<> mp;
    mp.addPatch(gsNurbsCreator<>::BSplineSquare());
    mp.computeTopology();

    gsParaview<real_t> pv;
    pv.options().setSwitch("singleFile", true);
    pv.options().setSwitch("elements",   true);
    pv.options().setSwitch("controlNet", true);
    pv.options().setSwitch("writePvd",   false);
    pv.write(mp, fn);

    CHECK(gsFileManager::fileExists(fn + ".vtu"));
    CHECK(gsFileManager::fileExists(fn + "_mesh.vtp"));
    CHECK(gsFileManager::fileExists(fn + "_cnet.vtp"));
    CHECK(!gsFileManager::fileExists(fn + ".pvd"));
}

TEST(SingleFileField_pointDataOnly_withExtras)
{
    // A gsField in singleFile mode must write SolutionField exactly once, as
    // PointData, with PatchID as cell data (one tuple per cell) -- and must
    // still produce the merged _mesh.vtp / _cnet.vtp like the gsMultiPatch
    // overload does.
    const std::string tmp = gsFileManager::getTempPath();
    if (tmp.empty()) return;

    const std::string fn = tmp + "gsParaview_field_sf_extras_test";

    // getTempPath() falls back to the working directory when TMPDIR is
    // unset, so outputs persist between runs: a negative check below would
    // otherwise report on a leftover file rather than on this run.
    std::remove((fn + ".pvd").c_str());
    std::remove((fn + ".vtu").c_str());
    std::remove((fn + "_mesh.vtp").c_str());
    std::remove((fn + "_cnet.vtp").c_str());

    gsMultiPatch<> mp;
    mp.addPatch(gsNurbsCreator<>::BSplineSquare());
    mp.addPatch(gsNurbsCreator<>::BSplineSquare(1, 1, 0));
    mp.computeTopology();

    gsField<> field(mp, mp);

    gsParaview<real_t> pv;
    pv.options().setSwitch("singleFile", true);
    pv.options().setSwitch("elements",   true);
    pv.options().setSwitch("controlNet", true);
    pv.write(field, fn);

    CHECK(gsFileManager::fileExists(fn + ".pvd"));
    CHECK(gsFileManager::fileExists(fn + ".vtu"));
    CHECK(gsFileManager::fileExists(fn + "_mesh.vtp"));
    CHECK(gsFileManager::fileExists(fn + "_cnet.vtp"));

    std::ifstream f(fn + ".vtu");
    std::string s((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());

    // SolutionField occurs exactly once as a named DataArray (as PointData);
    // a second, malformed CellData copy is the defect this test guards
    // against. (The attribute "Scalars=\"SolutionField\"" / "Vectors=..." on
    // the enclosing <PointData> tag also contains the substring, so counting
    // Name="SolutionField" is what distinguishes one array from a duplicate.)
    const std::string nameTag = "Name=\"SolutionField\"";
    std::size_t count = 0;
    std::size_t pos = 0;
    while ((pos = s.find(nameTag, pos)) != std::string::npos)
    {
        ++count;
        pos += nameTag.size();
    }
    CHECK(count > 0);
    CHECK_EQUAL(std::size_t(1), count);

    CHECK(s.find("<CellData Scalars=\"SolutionField\">") == std::string::npos);

    const std::string noc = "NumberOfCells=\"";
    std::size_t nocPos = s.find(noc);
    CHECK(nocPos != std::string::npos);
    std::size_t nStart = nocPos + noc.size();
    std::size_t nEnd = s.find('"', nStart);
    CHECK(nEnd != std::string::npos);
    const index_t nCells = atoi(s.substr(nStart, nEnd - nStart).c_str());
    CHECK(nCells > 0);

    const std::string patchIdTag = "Name=\"PatchID\"";
    std::size_t pidPos = s.find(patchIdTag);
    CHECK(pidPos != std::string::npos);
    std::size_t tagEnd = s.find('>', pidPos);
    CHECK(tagEnd != std::string::npos);
    std::size_t bodyStart = tagEnd + 1;
    std::size_t bodyEnd = s.find("</DataArray>", bodyStart);
    CHECK(bodyEnd != std::string::npos);
    const std::string body = s.substr(bodyStart, bodyEnd - bodyStart);

    std::istringstream iss(body);
    index_t tokenCount = 0;
    std::string token;
    while (iss >> token)
        ++tokenCount;

    CHECK(tokenCount > 0);
    CHECK_EQUAL(nCells, tokenCount);

    // The .pvd must actually reference the three files it was written
    // alongside, not merely leave them present on disk: mutating the
    // "skipPvd || extras" guard back to "skipPvd" at gsParaview.hpp:208
    // still writes _mesh.vtp / _cnet.vtp but stops listing them here.
    std::ifstream pvdFile(fn + ".pvd");
    std::string pvd((std::istreambuf_iterator<char>(pvdFile)), std::istreambuf_iterator<char>());

    const std::string dataSetTag = "<DataSet";
    std::size_t dsCount = 0;
    std::size_t dsPos = 0;
    while ((dsPos = pvd.find(dataSetTag, dsPos)) != std::string::npos)
    {
        ++dsCount;
        dsPos += dataSetTag.size();
    }
    CHECK(dsCount > 0);
    CHECK_EQUAL(std::size_t(3), dsCount);

    const std::string base = gsFileManager::getFilename(fn);
    CHECK(pvd.find(base + ".vtu") != std::string::npos);
    CHECK(pvd.find(base + "_mesh.vtp") != std::string::npos);
    CHECK(pvd.find(base + "_cnet.vtp") != std::string::npos);
}

TEST(SingleFileField_noPvd)
{
    // "writePvd"=false on the gsField overload: the merged _mesh.vtp /
    // _cnet.vtp must still be written, and no collection file created.
    const std::string tmp = gsFileManager::getTempPath();
    if (tmp.empty()) return;

    const std::string fn = tmp + "gsParaview_field_sf_nopvd_test";

    // A stale .pvd from an earlier run would make the negative check below
    // pass or fail for the wrong reason.
    std::remove((fn + ".pvd").c_str());
    std::remove((fn + ".vtu").c_str());
    std::remove((fn + "_mesh.vtp").c_str());
    std::remove((fn + "_cnet.vtp").c_str());

    gsMultiPatch<> mp;
    mp.addPatch(gsNurbsCreator<>::BSplineSquare());
    mp.addPatch(gsNurbsCreator<>::BSplineSquare(1, 1, 0));
    mp.computeTopology();

    gsField<> field(mp, mp);

    gsParaview<real_t> pv;
    pv.options().setSwitch("singleFile",  true);
    pv.options().setSwitch("elements",    true);
    pv.options().setSwitch("controlNet",  true);
    pv.options().setSwitch("writePvd",    false);
    pv.write(field, fn);

    CHECK(gsFileManager::fileExists(fn + ".vtu"));
    CHECK(gsFileManager::fileExists(fn + "_mesh.vtp"));
    CHECK(gsFileManager::fileExists(fn + "_cnet.vtp"));
    CHECK(!gsFileManager::fileExists(fn + ".pvd"));
}

TEST(TimeSteppingElementMesh_smoke)
{
    const std::string tmp = gsFileManager::getTempPath();
    if (tmp.empty()) return;

    const std::string fn = tmp + "gsParaview_ts_elements_test";

    gsMultiPatch<> mp;
    mp.addPatch(gsNurbsCreator<>::BSplineSquare());
    gsField<> field(mp, mp);

    gsParaviewCollection<real_t> collection(fn);
    collection.options().setInt("numPoints", 64);
    collection.options().setSwitch("elements", true);
    collection.newTimeStep(mp, 0.0);
    collection.addField(field, "identity");
    collection.saveTimeStep();
    collection.save();

    CHECK(gsFileManager::fileExists(fn + ".pvd"));
    CHECK(gsFileManager::fileExists(fn + "_pvd/" + gsFileManager::getBasename(fn) + "_t0.000000_mesh0.vtp"));
}

} // SUITE
