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

#include <algorithm> // std::min_element
#include <cstdio> // std::remove
#include <fstream>
#include <vector>
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
    CHECK(pv.options().getSwitch("multiblock"));
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
        pv.options().setSwitch("multiblock", false);
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
        pv.options().setSwitch("multiblock", false);
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
    // With "multiblock" disabled, the element mesh and the control net are
    // merged over all patches into one file each, so the output count does not
    // grow with nPatches(). Two patches must still give exactly one _mesh and
    // one _cnet.
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
    pv.options().setSwitch("multiblock", false);
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
    pv.options().setSwitch("multiblock", false);
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
    // A gsField in single-file mode must write SolutionField exactly once, as
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
    pv.options().setSwitch("multiblock", false);
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
    pv.options().setSwitch("multiblock",  false);
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

// Reads a whole file into a string (empty if it cannot be opened).
static std::string readFile(const std::string& path)
{
    std::ifstream f(path.c_str());
    return std::string((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

// Number of grid points (a+1)(b+1)(c+1) of the WholeExtent="0 a 0 b 0 c" of a .vts file.
static index_t extentPoints(const std::string& s)
{
    const std::string tag = "WholeExtent=\"";
    const size_t p = s.find(tag);
    if (p == std::string::npos) return -1;
    std::istringstream iss(s.substr(p + tag.size(), s.find('"', p + tag.size()) - p - tag.size()));
    index_t lo, hi, n = 1;
    for (int d = 0; d < 3; ++d)
    {
        iss >> lo >> hi;
        n *= hi - lo + 1;
    }
    return n;
}

// Values of the ascii DataArray named label; empty if there is none.
static std::vector<double> dataArrayValues(const std::string& s, const std::string& label)
{
    std::vector<double> v;
    const size_t p = s.find("Name=\"" + label + "\"");
    if (p == std::string::npos) return v;
    const size_t b = s.find('>', p) + 1;
    const size_t e = s.find("</DataArray>", b);
    std::istringstream iss(s.substr(b, e - b));
    double x;
    while (iss >> x) v.push_back(x);
    return v;
}

// Writes the field x+2y (or u+2v) on two unit-square patches through a
// collection and checks the pieces; lo[k], hi[k] are the expected value
// ranges on patch k.
static void checkExprFieldCollection(const bool isParam, const std::string& base,
                                     const std::string& label,
                                     const double lo[2], const double hi[2])
{
    const std::string tmp = gsFileManager::getTempPath();
    if (tmp.empty()) return;

    const std::string dir = tmp + "gsParaview_exprfield_test/";
    const std::string fn  = dir + base;
    const std::string sub = fn + "_pvd/";
    std::string piece[2];
    for (int k = 0; k < 2; ++k)
    {
        std::ostringstream os;
        os << sub << base << "_t0.000000_patch" << k << ".vts";
        piece[k] = os.str();
        std::remove(piece[k].c_str());
    }
    std::remove((fn + ".pvd").c_str());

    gsMultiPatch<> mp;
    mp.addPatch(gsNurbsCreator<>::BSplineSquare());
    mp.addPatch(gsNurbsCreator<>::BSplineSquare(1, 1, 0));
    mp.computeTopology();

    gsFunctionExpr<> f("x + 2*y", 2);
    gsField<> field(mp, f, isParam);

    gsParaviewCollection<real_t> collection(fn);
    collection.options().setInt("numPoints", 64);
    collection.newTimeStep(mp, 0.0);
    collection.addField(field, label);
    collection.saveTimeStep();
    collection.save();

    CHECK(gsFileManager::fileExists(fn + ".pvd"));
    for (int k = 0; k < 2; ++k)
    {
        CHECK(gsFileManager::fileExists(piece[k]));
        const std::string s = readFile(piece[k]);
        const std::vector<double> v = dataArrayValues(s, label);
        CHECK(!v.empty());
        CHECK_EQUAL(extentPoints(s), static_cast<index_t>(v.size()));
        if (!v.empty())
        {
            CHECK_CLOSE(lo[k], *std::min_element(v.begin(), v.end()), 1e-4);
            CHECK_CLOSE(hi[k], *std::max_element(v.begin(), v.end()), 1e-4);
        }
    }

    for (int k = 0; k < 2; ++k)
        std::remove(piece[k].c_str());
    std::remove(sub.c_str());
    std::remove((fn + ".pvd").c_str());
    std::remove(dir.c_str());
}

TEST(CollectionNonParametricExprField)
{
    // x+2y over the physical patches [0,1]x[0,1] and [1,2]x[0,1]
    const double lo[2] = {0.0, 1.0}, hi[2] = {3.0, 4.0};
    checkExprFieldCollection(false, "nonparam", "nonparam", lo, hi);
}

TEST(CollectionParametricExprField)
{
    // u+2v over the parameter domain [0,1]^2 of both patches
    const double lo[2] = {0.0, 0.0}, hi[2] = {3.0, 3.0};
    checkExprFieldCollection(true, "param", "param", lo, hi);
}

} // SUITE
