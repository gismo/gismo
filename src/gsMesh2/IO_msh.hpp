/** @file IO_msh.hpp

    @brief Gmsh (.msh) reader and writer for gsVolMesh.

    Reads the ASCII versions 2.2 and 4.1 of the format and writes ASCII 2.2,
    which every tool understands.  Binary files are rejected explicitly.

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
#include <map>
#include <sstream>
#include <string>

namespace gismo {

namespace internal {

/// advance \a in until the line that is exactly \a section, e.g. "$Nodes"
inline bool msh_seek(std::istream& in, const std::string& section)
{
    std::string line;
    while (std::getline(in, line))
    {
        while (!line.empty() && ('\r'==line[line.size()-1] || ' '==line[line.size()-1]))
            line.erase(line.size()-1);
        if (line == section) return true;
    }
    return false;
}

} // namespace internal

template <class Scalar>
bool read_msh(gsVolMesh<Scalar>& mesh, const std::string& filename)
{
    typedef typename gsVolMesh<Scalar>::Vertex Vertex;
    typedef typename gsVolMesh<Scalar>::Point  Point;

    std::ifstream in(filename.c_str());
    if (!in)
    {
        gsWarn << "read_msh: cannot open \"" << filename << "\".\n";
        return false;
    }

    // ---- header ---------------------------------------------------------
    if (!internal::msh_seek(in, "$MeshFormat"))
    {
        gsWarn << "read_msh: \"" << filename << "\" has no $MeshFormat section.\n";
        return false;
    }

    double version = 0; int filetype = 0, datasize = 0;
    in >> version >> filetype >> datasize;
    if (0 != filetype)
    {
        gsWarn << "read_msh: \"" << filename << "\" is a binary Gmsh file; only "
                  "ASCII is supported.\n";
        return false;
    }
    const bool v4 = (version >= 4.0);
    if (!v4 && version < 2.0)
    {
        gsWarn << "read_msh: unsupported Gmsh version " << version << ".\n";
        return false;
    }

    gsVolMesh<Scalar> out;
    std::map<long,Vertex> node;      // Gmsh node tags may be sparse

    // ---- nodes ----------------------------------------------------------
    if (!internal::msh_seek(in, "$Nodes"))
    {
        gsWarn << "read_msh: \"" << filename << "\" has no $Nodes section.\n";
        return false;
    }

    if (v4)
    {
        long nblocks=0, nnodes=0, mintag=0, maxtag=0;
        in >> nblocks >> nnodes >> mintag >> maxtag;
        out.reserve((unsigned)nnodes, 0, 0, 0);
        for (long b=0; b!=nblocks; ++b)
        {
            int dim=0, tag=0, param=0; long cnt=0;
            in >> dim >> tag >> param >> cnt;
            std::vector<long> tags((size_t)cnt);
            for (long i=0; i!=cnt; ++i) in >> tags[(size_t)i];
            for (long i=0; i!=cnt; ++i)
            {
                double x,y,z; in >> x >> y >> z;
                node[tags[(size_t)i]] =
                    out.add_vertex(Point((Scalar)x,(Scalar)y,(Scalar)z));
            }
        }
    }
    else
    {
        long nnodes=0;
        in >> nnodes;
        out.reserve((unsigned)nnodes, 0, 0, 0);
        for (long i=0; i!=nnodes; ++i)
        {
            long tag=0; double x,y,z;
            in >> tag >> x >> y >> z;
            node[tag] = out.add_vertex(Point((Scalar)x,(Scalar)y,(Scalar)z));
        }
    }

    // ---- elements -------------------------------------------------------
    if (!internal::msh_seek(in, "$Elements"))
    {
        gsWarn << "read_msh: \"" << filename << "\" has no $Elements section.\n";
        return false;
    }

    typename gsVolMesh<Scalar>::template Cell_property<index_t> ctag =
        out.template add_cell_property<index_t>("C:tag", 0);

    unsigned int skipped = 0, bad = 0;
    std::vector<Vertex> cv;

    // read one element given its Gmsh type and its node tags
    struct Local
    {
        static bool add(gsVolMesh<Scalar>& m, std::map<long,Vertex>& node,
                        int gtype, const std::vector<long>& ids,
                        std::vector<Vertex>& cv,
                        typename gsVolMesh<Scalar>::template Cell_property<index_t>& ctag,
                        index_t tag, unsigned int& skipped, unsigned int& bad)
        {
            const int vt = vtk_type_of_gmsh(gtype);
            if (0 == vt) { ++skipped; return true; }     // 0D/1D/2D marker

            cv.clear();
            for (size_t i=0; i!=ids.size(); ++i)
            {
                const typename std::map<long,Vertex>::const_iterator it = node.find(ids[i]);
                if (node.end() == it) { ++bad; return true; }
                cv.push_back(it->second);
            }
            const typename gsVolMesh<Scalar>::Cell c = add_typed_cell(m, vt, cv);
            if (!c.is_valid()) { ++bad; return true; }
            ctag[c] = tag;
            return true;
        }
    };

    if (v4)
    {
        long nblocks=0, nelem=0, mintag=0, maxtag=0;
        in >> nblocks >> nelem >> mintag >> maxtag;
        for (long b=0; b!=nblocks; ++b)
        {
            int dim=0, tag=0, gtype=0; long cnt=0;
            in >> dim >> tag >> gtype >> cnt;
            const int vt = vtk_type_of_gmsh(gtype);
            const unsigned int nn = (0!=vt) ? vtk_type_size(vt) : 0;

            for (long i=0; i!=cnt; ++i)
            {
                long etag=0; in >> etag;
                if (0 == vt)
                {
                    // an element we do not turn into a cell: skip its nodes
                    std::string rest; std::getline(in, rest);
                    ++skipped;
                    continue;
                }
                std::vector<long> ids(nn);
                for (unsigned int k=0; k!=nn; ++k) in >> ids[k];
                Local::add(out, node, gtype, ids, cv, ctag, (index_t)tag, skipped, bad);
            }
        }
    }
    else
    {
        long nelem=0;
        in >> nelem;
        for (long i=0; i!=nelem; ++i)
        {
            long etag=0; int gtype=0, ntags=0;
            in >> etag >> gtype >> ntags;
            index_t first = 0;
            for (int k=0; k!=ntags; ++k)
            {
                long t=0; in >> t;
                if (0==k) first = (index_t)t;
            }
            const int vt = vtk_type_of_gmsh(gtype);
            if (0 == vt)
            {
                std::string rest; std::getline(in, rest);
                ++skipped;
                continue;
            }
            const unsigned int nn = vtk_type_size(vt);
            std::vector<long> ids(nn);
            for (unsigned int k=0; k!=nn; ++k) in >> ids[k];
            Local::add(out, node, gtype, ids, cv, ctag, first, skipped, bad);
        }
    }

    if (0 == out.n_cells())
    {
        gsWarn << "read_msh: \"" << filename << "\" contains no 3D elements.\n";
        return false;
    }
    if (0 != bad)
        gsWarn << "read_msh: skipped " << bad << " malformed element(s).\n";
    if (0 != skipped)
        gsInfo << "read_msh: ignored " << skipped
               << " element(s) of dimension < 3 (boundary markers).\n";

    const unsigned int inv = internal::count_inverted(out);
    if (0 != inv)
        gsWarn << "read_msh: " << inv << " of " << out.n_cells()
               << " cells have non-positive volume.\n";

    mesh = give(out);
    return true;
}

template <class Scalar>
bool write_msh(const gsVolMesh<Scalar>& mesh, const std::string& filename)
{
    std::ofstream out(filename.c_str());
    if (!out)
    {
        gsWarn << "write_msh: cannot open \"" << filename << "\" for writing.\n";
        return false;
    }
    out.precision(16);

    out << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n";

    // Gmsh node tags are 1-based and the mesh may carry deleted vertices, so
    // the tags are assigned here rather than taken from the handles
    std::map<int,long> tag;
    long n = 0;
    out << "$Nodes\n" << mesh.n_vertices() << "\n";
    for (auto v : mesh.vertices())
    {
        tag[v.idx()] = ++n;
        const typename gsVolMesh<Scalar>::Point & p = mesh.position(v);
        out << n << " " << p[0] << " " << p[1] << " " << p[2] << "\n";
    }
    out << "$EndNodes\n";

    typename gsVolMesh<Scalar>::template Cell_property<index_t> ctag =
        mesh.template get_cell_property<index_t>("C:tag");

    std::vector<typename gsVolMesh<Scalar>::Vertex> ord;
    std::ostringstream body;
    long ne = 0, unsupported = 0;
    for (auto c : mesh.cells())
    {
        const int vt = mesh.cell_vtk_order(c, ord);
        const int gt = gmsh_type_of_vtk(vt);
        if (0 == gt) { ++unsupported; continue; }

        const index_t phys = ctag ? ctag[c] : 0;
        body << ++ne << " " << gt << " 2 " << phys << " " << phys;
        for (size_t i=0; i!=ord.size(); ++i) body << " " << tag[ord[i].idx()];
        body << "\n";
    }

    if (0 != unsupported)
    {
        gsWarn << "write_msh: " << unsupported << " cell(s) are not tetrahedra, "
                  "hexahedra, prisms or pyramids; the Gmsh format cannot express "
                  "them and they were not written.\n";
        if (0 == ne) return false;
    }

    out << "$Elements\n" << ne << "\n" << body.str() << "$EndElements\n";

    return out.good();
}

} // namespace gismo
