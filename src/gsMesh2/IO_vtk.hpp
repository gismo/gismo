/** @file IO_vtk.hpp

    @brief Legacy VTK (.vtk) reader and writer for gsVolMesh.

    Handles the ASCII DATASET UNSTRUCTURED_GRID flavour.  The legacy format has
    no unambiguous way of describing a general polyhedron, so only the four
    standard cell types are supported; anything else is reported rather than
    guessed at.  Use .vtu for polyhedral cells.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
*/

#pragma once

#include <gsMesh2/IO_vol.hpp>

#include <gsCore/gsDebug.h>

#include <fstream>
#include <sstream>
#include <string>

namespace gismo {

template <class Scalar>
bool read_vtk(gsVolMesh<Scalar>& mesh, const std::string& filename)
{
    typedef typename gsVolMesh<Scalar>::Vertex Vertex;
    typedef typename gsVolMesh<Scalar>::Point  Point;

    std::ifstream in(filename.c_str());
    if (!in)
    {
        gsWarn << "read_vtk: cannot open \"" << filename << "\".\n";
        return false;
    }

    std::string line;
    std::getline(in, line);                                   // # vtk DataFile
    if (std::string::npos == line.find("vtk"))
    {
        gsWarn << "read_vtk: \"" << filename << "\" is not a VTK file.\n";
        return false;
    }
    std::getline(in, line);                                   // title

    std::getline(in, line);                                   // ASCII | BINARY
    if (std::string::npos != line.find("BINARY"))
    {
        gsWarn << "read_vtk: \"" << filename << "\" is BINARY; only ASCII legacy "
                  "VTK is supported.\n";
        return false;
    }

    gsVolMesh<Scalar> out;
    std::vector<Vertex> pts;
    std::vector< std::vector<long> > conn;
    std::vector<int> types;

    std::string tok;
    while (in >> tok)
    {
        if ("DATASET" == tok)
        {
            std::string kind; in >> kind;
            if ("UNSTRUCTURED_GRID" != kind)
            {
                gsWarn << "read_vtk: dataset type " << kind
                       << " is not supported; expected UNSTRUCTURED_GRID.\n";
                return false;
            }
        }
        else if ("POINTS" == tok)
        {
            long n = 0; std::string dtype;
            in >> n >> dtype;
            pts.reserve((size_t)n);
            for (long i=0; i!=n; ++i)
            {
                double x,y,z; in >> x >> y >> z;
                pts.push_back(out.add_vertex(Point((Scalar)x,(Scalar)y,(Scalar)z)));
            }
        }
        else if ("CELLS" == tok)
        {
            long n=0, total=0;
            in >> n >> total;
            conn.resize((size_t)n);
            for (long i=0; i!=n; ++i)
            {
                long k=0; in >> k;
                conn[(size_t)i].resize((size_t)k);
                for (long j=0; j!=k; ++j) in >> conn[(size_t)i][(size_t)j];
            }
        }
        else if ("CELL_TYPES" == tok)
        {
            long n=0; in >> n;
            types.resize((size_t)n);
            for (long i=0; i!=n; ++i) in >> types[(size_t)i];
        }
        else if ("POINT_DATA" == tok || "CELL_DATA" == tok)
        {
            break;                       // field data is not read
        }
    }

    if (conn.size() != types.size())
    {
        gsWarn << "read_vtk: CELLS and CELL_TYPES disagree ("
               << conn.size() << " vs " << types.size() << ").\n";
        return false;
    }

    unsigned int bad = 0;
    std::vector<Vertex> cv;
    for (size_t i=0; i!=conn.size(); ++i)
    {
        if (0 == vtk_type_size(types[i])) { ++bad; continue; }

        cv.clear();
        bool ok = true;
        for (size_t j=0; j!=conn[i].size(); ++j)
        {
            const long id = conn[i][j];
            if (id < 0 || id >= (long)pts.size()) { ok = false; break; }
            cv.push_back(pts[(size_t)id]);
        }
        if (!ok || !add_typed_cell(out, types[i], cv).is_valid()) ++bad;
    }

    if (0 == out.n_cells())
    {
        gsWarn << "read_vtk: \"" << filename << "\" contains no supported 3D cells.\n";
        return false;
    }
    if (0 != bad)
        gsWarn << "read_vtk: skipped " << bad << " cell(s) of an unsupported type; "
                  "the legacy format cannot express polyhedra, use .vtu.\n";

    const unsigned int inv = internal::count_inverted(out);
    if (0 != inv)
        gsWarn << "read_vtk: " << inv << " of " << out.n_cells()
               << " cells have non-positive volume.\n";

    mesh = give(out);
    return true;
}

template <class Scalar>
bool write_vtk(const gsVolMesh<Scalar>& mesh, const std::string& filename)
{
    std::ofstream out(filename.c_str());
    if (!out)
    {
        gsWarn << "write_vtk: cannot open \"" << filename << "\" for writing.\n";
        return false;
    }
    out.precision(16);

    std::map<int,long> id;
    long n = 0;
    for (auto v : mesh.vertices()) id[v.idx()] = n++;

    std::vector<typename gsVolMesh<Scalar>::Vertex> ord;
    std::ostringstream cells, ctypes;
    long ncells = 0, total = 0, unsupported = 0;
    for (auto c : mesh.cells())
    {
        const int vt = mesh.cell_vtk_order(c, ord);
        if (0 == vtk_type_size(vt)) { ++unsupported; continue; }

        cells << ord.size();
        for (size_t i=0; i!=ord.size(); ++i) cells << " " << id[ord[i].idx()];
        cells << "\n";
        ctypes << vt << "\n";
        total += 1 + (long)ord.size();
        ++ncells;
    }

    if (0 != unsupported)
    {
        gsWarn << "write_vtk: " << unsupported << " polyhedral cell(s) cannot be "
                  "written to the legacy VTK format; use .vtu instead.\n";
        if (0 == ncells) return false;
    }

    out << "# vtk DataFile Version 3.0\n"
        << "gsVolMesh\nASCII\nDATASET UNSTRUCTURED_GRID\n"
        << "POINTS " << n << " double\n";
    for (auto v : mesh.vertices())
    {
        const typename gsVolMesh<Scalar>::Point & p = mesh.position(v);
        out << p[0] << " " << p[1] << " " << p[2] << "\n";
    }
    out << "CELLS " << ncells << " " << total << "\n" << cells.str()
        << "CELL_TYPES " << ncells << "\n" << ctypes.str();

    return out.good();
}

} // namespace gismo
