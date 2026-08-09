/** @file gsVolMesh.hpp

    @brief Implementation of gsVolMesh, the geometry layer of the half-face
    volume mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
*/

#pragma once

#include <gsMesh2/gsVolMesh.h>

#include <gsCore/gsDebug.h>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <map>
#include <sstream>

namespace gismo {

// ===========================================================================
//  construction, assignment
// ===========================================================================

template <class Scalar>
gsVolMesh<Scalar>::
gsVolMesh() : Base()
{
    // the topology properties are allocated by the base class constructor;
    // here we only add the geometry
    vpoint_ = add_vertex_property<Point>("V:point", Point(0,0,0));
}

template <class Scalar>
gsVolMesh<Scalar>::
gsVolMesh(const gsMatrix<Scalar> & pts) : gsVolMesh()
{
    GISMO_ASSERT(3 == pts.rows(),
                 "gsVolMesh: expecting three coordinates per point, got "<<pts.rows());
    for (index_t i = 0; i != pts.cols(); ++i)
        this->add_vertex(Point(pts(0,i), pts(1,i), pts(2,i)));
}

template <class Scalar>
gsVolMesh<Scalar>::
~gsVolMesh()
{
}

template <class Scalar>
gsVolMesh<Scalar>&
gsVolMesh<Scalar>::
operator=(const gsVolMesh<Scalar>& rhs)
{
    if (this != &rhs)
    {
        Base::operator=(rhs);
        // property handles contain pointers, have to be reassigned
        vpoint_ = get_vertex_property<Point>("V:point");
    }
    return *this;
}

template <class Scalar>
gsVolMesh<Scalar>&
gsVolMesh<Scalar>::
operator=(gsVolMesh<Scalar>&& rhs) noexcept
{
    if (this != &rhs)
    {
        Base::operator=(std::move(rhs));
        vpoint_ = get_vertex_property<Point>("V:point");
    }
    return *this;
}

template <class Scalar>
gsVolMesh<Scalar>&
gsVolMesh<Scalar>::
assign(const gsVolMesh<Scalar>& rhs)
{
    if (this != &rhs)
    {
        // clears the property containers and re-creates the topology properties
        Base::assign(rhs);

        // re-create and copy the geometry
        vpoint_ = add_vertex_property<Point>("V:point", Point(0,0,0));
        vpoint_.array() = rhs.vpoint_.array();
    }
    return *this;
}

template <class Scalar>
typename gsVolMesh<Scalar>::Vertex
gsVolMesh<Scalar>::
add_vertex(const Point& p)
{
    const Vertex v = Base::add_vertex();
    vpoint_[v] = p;
    return v;
}

// ===========================================================================
//  geometry
// ===========================================================================

template <class Scalar>
Scalar
gsVolMesh<Scalar>::
edge_length(Edge e) const
{
    return (vpoint_[vertex(e,0)] - vpoint_[vertex(e,1)]).norm();
}

template <class Scalar>
typename gsVolMesh<Scalar>::Point
gsVolMesh<Scalar>::
barycenter(Halfface hf) const
{
    Point c(0,0,0);
    unsigned int n = 0;
    for (auto v : vertices(hf)) { c += vpoint_[v]; ++n; }
    GISMO_ASSERT(0 != n, "gsVolMesh::barycenter: empty half-face");
    return c / (Scalar)n;
}

template <class Scalar>
typename gsVolMesh<Scalar>::Point
gsVolMesh<Scalar>::
barycenter(Cell c) const
{
    Point b(0,0,0);
    unsigned int n = 0;
    for (auto cn : corners(c)) { b += vpoint_[vertex(cn)]; ++n; }
    GISMO_ASSERT(0 != n, "gsVolMesh::barycenter: empty cell");
    return b / (Scalar)n;
}

template <class Scalar>
typename gsVolMesh<Scalar>::Normal
gsVolMesh<Scalar>::
normal(Halfface hf) const
{
    // Newell's method: well defined also for a non-planar polygon, and its
    // magnitude is twice the projected area
    Normal nrm(0,0,0);
    const Halfedge start = halfedge(hf);
    Halfedge h = start;
    do
    {
        const Point & a = vpoint_[from_vertex(h)];
        const Point & b = vpoint_[to_vertex(h)];
        nrm[0] += (a[1]-b[1]) * (a[2]+b[2]);
        nrm[1] += (a[2]-b[2]) * (a[0]+b[0]);
        nrm[2] += (a[0]-b[0]) * (a[1]+b[1]);
        h = next_halfedge(h);
    } while (h != start);

    return nrm;
}

template <class Scalar>
Scalar
gsVolMesh<Scalar>::
volume(Cell c) const
{
    // divergence theorem: 6V = sum over the outward boundary triangles of
    // a . (b x c); each half-face is triangulated as a fan from its first vertex
    Scalar six_v = (Scalar)0;

    for (auto hf : halffaces(c))
    {
        const Halfedge start = halfedge(hf);
        const Point & a = vpoint_[from_vertex(start)];

        Halfedge h = next_halfedge(start);
        while (next_halfedge(h) != start)
        {
            const Point & b = vpoint_[from_vertex(h)];
            const Point & d = vpoint_[to_vertex(h)];
            six_v += a.dot(b.cross(d));
            h = next_halfedge(h);
        }
    }

    return six_v / (Scalar)6;
}

template <class Scalar>
Scalar
gsVolMesh<Scalar>::
volume() const
{
    Scalar v = (Scalar)0;
    for (auto c : cells()) v += volume(c);
    return v;
}

// ===========================================================================
//  surface views
// ===========================================================================

template <class Scalar>
typename gsVolMesh<Scalar>::Surface
gsVolMesh<Scalar>::
boundary_mesh() const
{
    Surface out;
    typename Surface::template Vertex_property<index_t> back =
        out.template add_vertex_property<index_t>("v:volvertex", -1);

    std::map<int, typename Surface::Vertex> known;
    std::vector<typename Surface::Vertex> loop;

    for (auto f : faces())
    {
        if (!is_boundary(f)) continue;

        loop.clear();
        for (auto v : vertices(halfface(f,0)))
        {
            typename std::map<int, typename Surface::Vertex>::iterator
                it = known.find(v.idx());
            if (known.end() == it)
            {
                const typename Surface::Vertex sv = out.add_vertex(vpoint_[v]);
                back[sv] = v.idx();
                it = known.insert(std::make_pair(v.idx(), sv)).first;
            }
            loop.push_back(it->second);
        }
        out.add_face(loop);
    }

    return out;
}

template <class Scalar>
typename gsVolMesh<Scalar>::Surface
gsVolMesh<Scalar>::
cell_mesh(Cell c) const
{
    Surface out;
    typename Surface::template Vertex_property<index_t> back =
        out.template add_vertex_property<index_t>("v:volvertex", -1);

    std::map<int, typename Surface::Vertex> known;
    for (auto cn : corners(c))
    {
        const Vertex v = vertex(cn);
        const typename Surface::Vertex sv = out.add_vertex(vpoint_[v]);
        back[sv] = v.idx();
        known[v.idx()] = sv;
    }

    std::vector<typename Surface::Vertex> loop;
    for (auto hf : halffaces(c))
    {
        loop.clear();
        for (auto v : vertices(hf)) loop.push_back(known[v.idx()]);
        out.add_face(loop);
    }

    return out;
}

// ===========================================================================
//  input / output
// ===========================================================================

template <class Scalar>
bool
gsVolMesh<Scalar>::
write(const std::string& filename) const
{
    const std::string::size_type dot = filename.rfind('.');
    if (std::string::npos == dot)
    {
        gsWarn << "gsVolMesh::write: no file extension in \""<<filename<<"\".\n";
        return false;
    }

    std::string ext = filename.substr(dot+1);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    if ("vtu" == ext) return write_vtu(filename);

    // every other format we know describes a surface, so write the boundary
    return boundary_mesh().write(filename);
}

template <class Scalar>
bool
gsVolMesh<Scalar>::
write_vtu(const std::string& filename) const
{
    std::ofstream out(filename.c_str());
    if (!out) return false;

    // VTK needs contiguous point ids, which the mesh only has after a garbage
    // collection; map them explicitly instead of requiring one
    std::map<int,index_t> vid;
    index_t nv = 0;
    for (auto v : vertices()) vid[v.idx()] = nv++;

    std::vector< std::vector<index_t> > conn;
    std::vector<int> types;

    for (auto c : cells())
    {
        std::vector<index_t> ids;
        for (auto cn : corners(c)) ids.push_back(vid[vertex(cn).idx()]);

        // VTK_TETRA / VTK_PYRAMID / VTK_WEDGE / VTK_HEXAHEDRON expect a fixed
        // corner ordering that the corner ring does not reproduce, so only the
        // general polyhedron cell type is emitted here
        conn.push_back(ids);
        types.push_back(42); // VTK_POLYHEDRON
    }

    out << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
        << "  <UnstructuredGrid>\n"
        << "    <Piece NumberOfPoints=\""<<nv<<"\" NumberOfCells=\""<<conn.size()<<"\">\n"
        << "      <Points>\n"
        << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (auto v : vertices())
        out << "          " << vpoint_[v][0] << " " << vpoint_[v][1] << " "
            << vpoint_[v][2] << "\n";
    out << "        </DataArray>\n"
        << "      </Points>\n"
        << "      <Cells>\n"
        << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (size_t i=0; i!=conn.size(); ++i)
    {
        out << "         ";
        for (size_t j=0; j!=conn[i].size(); ++j) out << " " << conn[i][j];
        out << "\n";
    }
    out << "        </DataArray>\n"
        << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n         ";
    {
        index_t off = 0;
        for (size_t i=0; i!=conn.size(); ++i) { off += (index_t)conn[i].size(); out << " " << off; }
    }
    out << "\n        </DataArray>\n"
        << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n         ";
    for (size_t i=0; i!=types.size(); ++i) out << " " << types[i];
    out << "\n        </DataArray>\n";

    // VTK_POLYHEDRON additionally needs the face stream of every cell
    out << "        <DataArray type=\"Int64\" Name=\"faces\" format=\"ascii\">\n";
    std::vector<index_t> faceoff;
    index_t stream = 0;
    for (auto c : cells())
    {
        std::ostringstream line;
        index_t nf = 0;
        for (auto hf : halffaces(c))
        {
            index_t nvf = 0;
            std::ostringstream lf;
            for (auto v : vertices(hf)) { lf << " " << vid[v.idx()]; ++nvf; }
            line << " " << nvf << lf.str();
            stream += 1 + nvf;
            ++nf;
        }
        out << "          " << nf << line.str() << "\n";
        stream += 1;
        faceoff.push_back(stream);
    }
    out << "        </DataArray>\n"
        << "        <DataArray type=\"Int64\" Name=\"faceoffsets\" format=\"ascii\">\n         ";
    for (size_t i=0; i!=faceoff.size(); ++i) out << " " << faceoff[i];
    out << "\n        </DataArray>\n"
        << "      </Cells>\n"
        << "    </Piece>\n"
        << "  </UnstructuredGrid>\n"
        << "</VTKFile>\n";

    return out.good();
}

} // namespace gismo
