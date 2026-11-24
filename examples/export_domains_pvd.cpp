/** @file export_domains_pvd.cpp

    @brief Utility to export G+Smo geometry XML files (single gsGeometry or
           gsMultiPatch, polynomial B-spline or rational NURBS) to ParaView
           .pvd collections with optional element mesh, control net, and
           geometric fields (measure, surface normal).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gismo.h>

#include <algorithm>
#include <dirent.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <vector>

using namespace gismo;

namespace {

void mkdirRecursive(const std::string& path)
{
    std::string cur;
    for (size_t i = 0; i < path.size(); ++i)
    {
        cur.push_back(path[i]);
        if (path[i] == '/' && cur.size() > 1)
            ::mkdir(cur.c_str(), 0755);
    }
    if (!cur.empty())
        ::mkdir(cur.c_str(), 0755);
}

std::string stripPathAndExt(const std::string& path)
{
    std::string base = path;
    const size_t slash = base.find_last_of("/\\");
    if (slash != std::string::npos)
        base = base.substr(slash + 1);
    const size_t dot = base.find_last_of('.');
    if (dot != std::string::npos)
        base = base.substr(0, dot);
    return base;
}

bool endsWith(const std::string& s, const std::string& suffix)
{
    return s.size() >= suffix.size() &&
           s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::vector<std::string> listXmlFiles(const std::string& dirPath)
{
    std::vector<std::string> files;
    DIR* dir = opendir(dirPath.c_str());
    if (!dir)
        return files;
    while (struct dirent* ent = readdir(dir))
    {
        std::string name(ent->d_name);
        if (name == "." || name == "..")
            continue;
        if (endsWith(name, ".xml"))
        {
            std::string full = dirPath;
            if (!full.empty() && full.back() != '/')
                full += '/';
            full += name;
            files.push_back(full);
        }
    }
    closedir(dir);
    std::sort(files.begin(), files.end());
    return files;
}

bool readAsMultiPatch(const std::string& filename, gsMultiPatch<>& mp)
{
    mp.clear();
    gsFileData<> fd(filename);
    if (fd.has<gsMultiPatch<>>())
    {
        mp = *fd.getFirst<gsMultiPatch<>>();
        return mp.nPatches() > 0;
    }
    if (fd.has<gsGeometry<>>())
    {
        std::vector<gsGeometry<>::uPtr> geoms = fd.getAll<gsGeometry<>>();
        for (size_t i = 0; i < geoms.size(); ++i)
            mp.addPatch(give(*geoms[i]));
        mp.computeTopology();
        return mp.nPatches() > 0;
    }
    return false;
}

bool exportOne(const std::string& xmlPath,
               const std::string& outDir,
               index_t numPoints,
               bool plotMesh,
               bool plotCtrlNet,
               bool withFields)
{
    gsMultiPatch<> mp;
    try
    {
        if (!readAsMultiPatch(xmlPath, mp))
        {
            gsWarn << "  [SKIP] No geometry found in: " << xmlPath << "\n";
            return false;
        }
    }
    catch (const std::exception& e)
    {
        gsWarn << "  [FAIL] Could not read " << xmlPath << ": " << e.what() << "\n";
        return false;
    }

    const std::string base = stripPathAndExt(xmlPath);
    std::string dir = outDir;
    if (!dir.empty() && dir.back() != '/')
        dir += '/';
    mkdirRecursive(dir);
    if (withFields)
        mkdirRecursive(dir + "_vtk");

    const std::string pvdBase = dir + base;

    const index_t nP = mp.nPatches();
    const short_t parDim = mp.parDim();
    const short_t geoDim = mp.targetDim();
    const bool isNurbs = (dynamic_cast<const gsTensorNurbs<2, real_t>*>(&mp.patch(0)) != nullptr);
    const short_t degU = mp.patch(0).basis().degree(0);
    const short_t degV = (parDim > 1) ? mp.patch(0).basis().degree(1) : 0;

    bool exported = false;
    if (withFields && parDim == 2)
    {
        try
        {
            gsMultiBasis<> mb(mp);
            gsExprEvaluator<> ev;
            ev.setIntegrationDomain(mb.domain());
            gsExprEvaluator<>::geometryMap G = ev.getMap(mp);

            gsParaviewCollection col(pvdBase + ".pvd", &ev);
            col.options().setInt("numPoints", numPoints);
            col.options().setInt("precision", 6);
            col.options().setSwitch("plotElements", plotMesh);
            col.options().setInt("plotElements.resolution", -1);
            col.options().setSwitch("plotControlNet", plotCtrlNet);
            col.options().setSwitch("makeSubfolder", true);
            col.options().setString("subfolder", "_vtk");

            col.newTimeStep(&mp);
            col.addField(meas(G), "Measure");
            if (geoDim == 3)
                col.addField(usn(G), "SurfaceNormal");
            col.saveTimeStep();
            col.save();
            exported = true;
        }
        catch (const std::exception& e)
        {
            gsWarn << "  [WARN] Field export failed for " << base
                   << " (" << e.what() << "); falling back to plain geometry export.\n";
        }
    }

    if (!exported)
    {
        try
        {
            gsWriteParaview(mp, pvdBase, static_cast<unsigned>(numPoints), plotMesh, plotCtrlNet);
            exported = true;
        }
        catch (const std::exception& e)
        {
            gsWarn << "  [FAIL] Export failed for " << base << ": " << e.what() << "\n";
            return false;
        }
    }

    gsInfo << "  [OK] " << std::left << std::setw(42) << (base + ".pvd")
           << " patches=" << std::setw(3) << nP
           << " dim=" << parDim << "->" << geoDim
           << " deg=(" << degU << "," << degV << ")"
           << " type=" << (isNurbs ? "NURBS" : "BSpline")
           << " ifaces=" << mp.nInterfaces() << "\n";
    return true;
}

} // anonymous namespace

int main(int argc, char* argv[])
{
    std::string inputFile("");
    std::string inputDir("");
    std::string outDir("paraview_domains");
    index_t numPoints = 1000;
    bool noMesh = false;
    bool ctrlNet = false;
    bool noFields = false;

    gsCmdLine cmd("Export G+Smo geometry XML file(s) to ParaView (.pvd).");
    cmd.addString("f", "file", "Single input geometry XML file.", inputFile);
    cmd.addString("d", "dir", "Directory of input geometry XML files to batch-export.", inputDir);
    cmd.addString("o", "outdir", "Output directory for .pvd files.", outDir);
    cmd.addInt("n", "points", "Sampling points per patch (default: 1000).", numPoints);
    cmd.addSwitch("nomesh", "Disable element/patch mesh lines.", noMesh);
    cmd.addSwitch("ctrlnet", "Include control net in output.", ctrlNet);
    cmd.addSwitch("nofields", "Export plain geometry without Measure/SurfaceNormal fields.", noFields);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    std::vector<std::string> targets;
    if (!inputFile.empty())
        targets.push_back(inputFile);
    if (!inputDir.empty())
    {
        std::vector<std::string> dirFiles = listXmlFiles(inputDir);
        targets.insert(targets.end(), dirFiles.begin(), dirFiles.end());
    }

    if (targets.empty())
    {
        gsInfo << "Usage:\n"
               << "  export_domains_pvd -f <file.xml> -o <outdir>\n"
               << "  export_domains_pvd -d <xml_dir>  -o <outdir>\n";
        return EXIT_FAILURE;
    }

    gsInfo << "Exporting " << targets.size() << " geometry file(s) to " << outDir << "/\n";
    index_t ok = 0;
    for (size_t i = 0; i < targets.size(); ++i)
    {
        if (exportOne(targets[i], outDir, numPoints, !noMesh, ctrlNet, !noFields))
            ++ok;
    }

    gsInfo << "Done: " << ok << " / " << targets.size() << " exported to " << outDir << "/\n";
    return (ok > 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
