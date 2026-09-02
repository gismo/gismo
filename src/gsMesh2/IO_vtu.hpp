/** @file IO_vtu.hpp

    @brief VTK XML unstructured grid (.vtu) reader and writer for gsVolMesh.

    Only inline ASCII DataArrays are handled; binary, appended and compressed
    payloads are rejected with a message rather than silently yielding an empty
    mesh.  Unlike the legacy format, .vtu can express a general polyhedron
    through the faces/faceoffsets streams, so every gsVolMesh survives a round
    trip here.

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

namespace internal {

/// value of attribute \a name inside the tag text \a tag, or an empty string
inline std::string xml_attr(const std::string& tag, const std::string& name)
{
    const std::string key = name + "=\"";
    const std::string::size_type a = tag.find(key);
    if (std::string::npos == a) return std::string();
    const std::string::size_type b = tag.find('"', a+key.size());
    if (std::string::npos == b) return std::string();
    return tag.substr(a+key.size(), b-a-key.size());
}

/** the inline text of the `<DataArray>` whose Name attribute is \a name.

    \param doc  the whole file
    \param from search from this offset, so that the Points and Cells sections
                can be scanned separately
    \param ok   set to false if the array exists but is not inline ASCII
*/
inline std::string xml_dataarray(const std::string& doc, const std::string& name,
                                 std::string::size_type from, bool& ok)
{
    std::string::size_type p = from;
    while (true)
    {
        p = doc.find("<DataArray", p);
        if (std::string::npos == p) return std::string();

        const std::string::size_type close = doc.find('>', p);
        if (std::string::npos == close) return std::string();

        const std::string tag = doc.substr(p, close-p);
        if (xml_attr(tag, "Name") == name)
        {
            const std::string fmt = xml_attr(tag, "format");
            if (!fmt.empty() && "ascii" != fmt)
            {
                gsWarn << "read_vtu: the DataArray \"" << name << "\" is stored as "
                       << fmt << "; only inline ascii is supported.\n";
                ok = false;
                return std::string();
            }
            const std::string::size_type end = doc.find("</DataArray>", close);
            if (std::string::npos == end) return std::string();
            return doc.substr(close+1, end-close-1);
        }
        p = close;
    }
}

} // namespace internal

template <class Scalar>
bool read_vtu(gsVolMesh<Scalar>& mesh, const std::string& filename)
{
    typedef typename gsVolMesh<Scalar>::Vertex Vertex;
    typedef typename gsVolMesh<Scalar>::Point  Point;

    std::ifstream in(filename.c_str());
    if (!in)
    {
        gsWarn << "read_vtu: cannot open \"" << filename << "\".\n";
        return false;
    }
    std::ostringstream buf; buf << in.rdbuf();
    const std::string doc = buf.str();

    if (std::string::npos != doc.find("<AppendedData"))
    {
        gsWarn << "read_vtu: \"" << filename << "\" uses appended data, which is "
                  "not supported; write it with format=\"ascii\".\n";
        return false;
    }
    if (std::string::npos == doc.find("UnstructuredGrid"))
    {
        gsWarn << "read_vtu: \"" << filename << "\" is not an UnstructuredGrid.\n";
        return false;
    }

    const std::string::size_type piece = doc.find("<Piece");
    if (std::string::npos == piece)
    {
        gsWarn << "read_vtu: \"" << filename << "\" has no <Piece>.\n";
        return false;
    }
    const std::string::size_type pts_at   = doc.find("<Points", piece);
    const std::string::size_type cells_at = doc.find("<Cells", piece);
    if (std::string::npos == pts_at || std::string::npos == cells_at)
    {
        gsWarn << "read_vtu: \"" << filename << "\" has no <Points> or <Cells>.\n";
        return false;
    }

    bool ok = true;
    // the Points DataArray usually has no Name, so it is read positionally
    std::string coords;
    {
        const std::string::size_type d = doc.find("<DataArray", pts_at);
        const std::string::size_type close = doc.find('>', d);
        const std::string::size_type end   = doc.find("</DataArray>", close);
        if (std::string::npos == d || std::string::npos == end)
        {
            gsWarn << "read_vtu: malformed <Points> section.\n";
            return false;
        }
        const std::string tag = doc.substr(d, close-d);
        const std::string fmt = internal::xml_attr(tag, "format");
        if (!fmt.empty() && "ascii" != fmt)
        {
            gsWarn << "read_vtu: the points are stored as " << fmt
                   << "; only inline ascii is supported.\n";
            return false;
        }
        coords = doc.substr(close+1, end-close-1);
    }

    const std::string sconn  = internal::xml_dataarray(doc, "connectivity", cells_at, ok);
    const std::string soff   = internal::xml_dataarray(doc, "offsets",      cells_at, ok);
    const std::string stypes = internal::xml_dataarray(doc, "types",        cells_at, ok);
    const std::string sfaces = internal::xml_dataarray(doc, "faces",        cells_at, ok);
    const std::string sfoff  = internal::xml_dataarray(doc, "faceoffsets",  cells_at, ok);
    if (!ok) return false;

    if (sconn.empty() || soff.empty() || stypes.empty())
    {
        gsWarn << "read_vtu: \"" << filename
               << "\" lacks connectivity, offsets or types.\n";
        return false;
    }

    gsVolMesh<Scalar> out;
    std::vector<Vertex> pts;
    {
        std::istringstream s(coords);
        double x,y,z;
        while (s >> x >> y >> z)
            pts.push_back(out.add_vertex(Point((Scalar)x,(Scalar)y,(Scalar)z)));
    }

    std::vector<long> conn, off, faces, foff;
    std::vector<int>  types;
    { std::istringstream s(sconn);  long v; while (s >> v) conn.push_back(v); }
    { std::istringstream s(soff);   long v; while (s >> v) off.push_back(v); }
    { std::istringstream s(stypes); int  v; while (s >> v) types.push_back(v); }
    if (!sfaces.empty()) { std::istringstream s(sfaces); long v; while (s >> v) faces.push_back(v); }
    if (!sfoff.empty())  { std::istringstream s(sfoff);  long v; while (s >> v) foff.push_back(v); }

    if (off.size() != types.size())
    {
        gsWarn << "read_vtu: offsets and types disagree ("
               << off.size() << " vs " << types.size() << ").\n";
        return false;
    }

    unsigned int bad = 0;
    std::vector<Vertex> cv;
    std::vector< std::vector<Vertex> > loops;

    for (size_t i=0; i!=types.size(); ++i)
    {
        const long begin = (0==i) ? 0 : off[i-1];
        const long end   = off[i];
        if (begin < 0 || end > (long)conn.size() || end < begin) { ++bad; continue; }

        if (42 == types[i])                                  // VTK_POLYHEDRON
        {
            if (foff.empty() || i >= foff.size()) { ++bad; continue; }
            const long fbegin = (0==i) ? 0 : (foff[i-1] < 0 ? 0 : foff[i-1]);
            // faceoffsets is -1 for non-polyhedral cells, so walk back for the
            // last real offset before this one
            long start = 0;
            for (long k=(long)i-1; k>=0; --k)
                if (foff[(size_t)k] >= 0) { start = foff[(size_t)k]; break; }
            GISMO_UNUSED(fbegin);

            if (start < 0 || start >= (long)faces.size()) { ++bad; continue; }
            long at = start;
            const long nf = faces[at++];
            loops.clear();
            bool okc = true;
            for (long f=0; f!=nf && okc; ++f)
            {
                if (at >= (long)faces.size()) { okc=false; break; }
                const long nvf = faces[at++];
                std::vector<Vertex> loop;
                for (long k=0; k!=nvf; ++k)
                {
                    if (at >= (long)faces.size()) { okc=false; break; }
                    const long id = faces[at++];
                    if (id < 0 || id >= (long)pts.size()) { okc=false; break; }
                    loop.push_back(pts[(size_t)id]);
                }
                loops.push_back(loop);
            }
            if (!okc || !out.add_cell(loops).is_valid()) ++bad;
            continue;
        }

        if (0 == vtk_type_size(types[i])) { ++bad; continue; }

        cv.clear();
        bool okc = true;
        for (long k=begin; k!=end; ++k)
        {
            const long id = conn[(size_t)k];
            if (id < 0 || id >= (long)pts.size()) { okc = false; break; }
            cv.push_back(pts[(size_t)id]);
        }
        if (!okc || !add_typed_cell(out, types[i], cv).is_valid()) ++bad;
    }

    if (0 == out.n_cells())
    {
        gsWarn << "read_vtu: \"" << filename << "\" contains no readable 3D cells.\n";
        return false;
    }
    if (0 != bad)
        gsWarn << "read_vtu: skipped " << bad << " cell(s) that could not be read.\n";

    const unsigned int inv = internal::count_inverted(out);
    if (0 != inv)
        gsWarn << "read_vtu: " << inv << " of " << out.n_cells()
               << " cells have non-positive volume.\n";

    mesh = give(out);
    return true;
}

template <class Scalar>
bool write_vtu(const gsVolMesh<Scalar>& mesh, const std::string& filename)
{
    std::ofstream out(filename.c_str());
    if (!out)
    {
        gsWarn << "write_vtu: cannot open \"" << filename << "\" for writing.\n";
        return false;
    }
    out.precision(16);

    // VTK needs contiguous point ids, which the mesh only has after a garbage
    // collection; map them explicitly rather than requiring one
    std::map<int,long> id;
    long nv = 0;
    for (auto v : mesh.vertices()) id[v.idx()] = nv++;

    std::vector<typename gsVolMesh<Scalar>::Vertex> ord;
    std::ostringstream sconn, soff, stypes, sfaces, sfoff;
    long ncells = 0, offset = 0, faceoffset = 0;
    bool any_poly = false;

    for (auto c : mesh.cells())
    {
        const int vt = mesh.cell_vtk_order(c, ord);

        sconn << "         ";
        for (size_t i=0; i!=ord.size(); ++i) sconn << " " << id[ord[i].idx()];
        sconn << "\n";
        offset += (long)ord.size();
        soff << " " << offset;
        stypes << " " << vt;
        ++ncells;

        if (42 != vt) { sfoff << " -1"; continue; }

        // a polyhedron additionally needs its face stream
        any_poly = true;
        std::ostringstream one;
        long nf = 0, nentries = 0;
        for (auto hf : mesh.halffaces(c))
        {
            long nvf = 0;
            std::ostringstream lf;
            for (auto v : mesh.vertices(hf)) { lf << " " << id[v.idx()]; ++nvf; }
            one << " " << nvf << lf.str();
            nentries += 1 + nvf;
            ++nf;
        }
        sfaces << "          " << nf << one.str() << "\n";
        faceoffset += 1 + nentries;
        sfoff << " " << faceoffset;
    }

    out << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
        << "  <UnstructuredGrid>\n"
        << "    <Piece NumberOfPoints=\""<<nv<<"\" NumberOfCells=\""<<ncells<<"\">\n"
        << "      <Points>\n"
        << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (auto v : mesh.vertices())
    {
        const typename gsVolMesh<Scalar>::Point & p = mesh.position(v);
        out << "          " << p[0] << " " << p[1] << " " << p[2] << "\n";
    }
    out << "        </DataArray>\n      </Points>\n      <Cells>\n"
        << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n"
        << sconn.str()
        << "        </DataArray>\n"
        << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n         "
        << soff.str() << "\n        </DataArray>\n"
        << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n         "
        << stypes.str() << "\n        </DataArray>\n";

    if (any_poly)
        out << "        <DataArray type=\"Int64\" Name=\"faces\" format=\"ascii\">\n"
            << sfaces.str()
            << "        </DataArray>\n"
            << "        <DataArray type=\"Int64\" Name=\"faceoffsets\" format=\"ascii\">\n         "
            << sfoff.str() << "\n        </DataArray>\n";

    out << "      </Cells>\n    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n";

    return out.good();
}

} // namespace gismo
