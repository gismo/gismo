/** @file gsVolCatmullClark.hpp

    @brief Implementation of Catmull-Clark subdivision on a volumetric mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
*/

#pragma once

#include <gsMesh2/gsSubdivisionSchemes/gsVolCatmullClark.h>

#include <gsCore/gsDebug.h>

#include <algorithm>
#include <vector>

namespace gismo {

template <class Scalar>
typename gsVolSubdivisionScheme<Scalar>::gsSubdivisionMeshValidity
gsVolCatmullClark<Scalar>::check_mesh()
{
    if (nullptr == this->m_mesh || 0 == this->m_mesh->n_cells())
        return Base::INVALID;
    return Base::VALID;
}

template <class Scalar>
void
gsVolCatmullClark<Scalar>::
subdivide_impl()
{
    gsVolMesh<Scalar> & m = this->mesh();
    if (0 == m.n_cells()) return;

    gsVolMesh<Scalar> out;
    out.reserve(m.n_vertices() + m.n_edges() + m.n_faces() + m.n_cells(),
                0, 0, m.n_corners());

    // the four point maps: properties on the old mesh holding handles into the
    // new one.  Both meshes use the same handle type, so no conversion is
    // needed; the maps disappear with the old mesh at the end.
    typename gsVolMesh<Scalar>::template Cell_property<Vertex> cellpt =
        m.template cell_property<Vertex>("C:cc-point", Vertex());
    typename gsVolMesh<Scalar>::template Face_property<Vertex> facept =
        m.template face_property<Vertex>("F:cc-point", Vertex());
    typename gsVolMesh<Scalar>::template Edge_property<Vertex> edgept =
        m.template edge_property<Vertex>("E:cc-point", Vertex());
    typename gsVolMesh<Scalar>::template Vertex_property<Vertex> vertpt =
        m.template vertex_property<Vertex>("V:cc-point", Vertex());

    // ---- cell points: the centroid of the cell's vertices ----------------
    for (auto c : m.cells())
        cellpt[c] = out.add_vertex(m.barycenter(c));

    // ---- face points: the face's vertices together with the cell points of
    //      the cells sharing it ------------------------------------------
    for (auto f : m.faces())
    {
        Point p(0,0,0);
        unsigned int n = 0;
        for (auto v : m.vertices(m.halfface(f,0))) { p += m.position(v); ++n; }
        for (index_t i = 0; i != 2; ++i)
        {
            const Cell c = m.cell(f,i);
            if (c.is_valid()) { p += out.position(cellpt[c]); ++n; }
        }
        facept[f] = out.add_vertex(p / (Scalar)n);
    }

    // ---- edge points: the two endpoints, the face points of the incident
    //      faces and the cell points of the incident cells ----------------
    std::vector<Face> efaces;
    for (auto e : m.edges())
    {
        Point p = m.position(m.vertex(e,0)) + m.position(m.vertex(e,1));
        unsigned int n = 2;

        // the radial orbit reports one dart per incident cell; each of those
        // cells contributes two faces at the edge, so collect and de-duplicate
        efaces.clear();
        for (auto h : m.halfedges(e))
        {
            efaces.push_back(m.face(m.halfface(h)));
            efaces.push_back(m.face(m.halfface(m.opposite_halfedge(h))));
            p += out.position(cellpt[m.cell(m.halfface(h))]);
            ++n;
        }
        std::sort(efaces.begin(), efaces.end());
        efaces.erase(std::unique(efaces.begin(), efaces.end()), efaces.end());

        for (size_t i = 0; i != efaces.size(); ++i)
        { p += out.position(facept[efaces[i]]); ++n; }

        edgept[e] = out.add_vertex(p / (Scalar)n);
    }

    // ---- vertex points: (A + 3B + 3C + V)/8, with A the average of the
    //      incident cell points, B of the incident face points and C of the
    //      incident edge midpoints ---------------------------------------
    for (auto v : m.vertices())
    {
        if (m.is_isolated(v)) continue;

        Point A(0,0,0); unsigned int na = 0;
        for (auto c : m.cells(v)) { A += out.position(cellpt[c]); ++na; }

        Point B(0,0,0); unsigned int nb = 0;
        const std::vector<Face> vf = m.faces(v);
        for (size_t i = 0; i != vf.size(); ++i)
        { B += out.position(facept[vf[i]]); ++nb; }

        Point C(0,0,0); unsigned int nc = 0;
        const std::vector<Edge> ve = m.edges(v);
        for (size_t i = 0; i != ve.size(); ++i)
        {
            C += (Scalar)0.5 * (m.position(m.vertex(ve[i],0)) +
                                m.position(m.vertex(ve[i],1)));
            ++nc;
        }

        GISMO_ASSERT(0!=na && 0!=nb && 0!=nc,
                     "gsVolCatmullClark: vertex "<<v<<" has an empty star");

        const Point p = ( A/(Scalar)na
                        + (Scalar)3 * B/(Scalar)nb
                        + (Scalar)3 * C/(Scalar)nc
                        + m.position(v) ) / (Scalar)8;
        vertpt[v] = out.add_vertex(p);
    }

    // ---- one new cell per (cell, corner) pair ---------------------------
    //
    // Walking the darts around a corner inside its cell gives the cyclic
    // sequence of edges e_i and faces f_i at that corner, with f_i bounded by
    // e_i and e_{i+1}.  The new cell is bounded by k quads around the vertex
    // and k quads around the cell point; every directed edge of that list
    // occurs once in each direction, which is what add_cell() verifies.
    std::vector<Edge> ce;
    std::vector<Face> cf;
    std::vector< std::vector<Vertex> > faces;

    for (auto c : m.cells())
    {
        const Vertex C = cellpt[c];

        for (auto cn : m.corners(c))
        {
            const Vertex v = m.vertex(cn);

            ce.clear();
            cf.clear();
            for (auto h : m.halfedges(cn))
            {
                ce.push_back(m.edge(h));
                cf.push_back(m.face(m.halfface(h)));
            }
            const size_t k = ce.size();
            GISMO_ASSERT(k >= 3, "gsVolCatmullClark: corner of valence "<<k);

            faces.clear();
            for (size_t i = 0; i != k; ++i)         // around the vertex
            {
                std::vector<Vertex> q;
                q.push_back(vertpt[v]);
                q.push_back(edgept[ce[i]]);
                q.push_back(facept[cf[i]]);
                q.push_back(edgept[ce[(i+1)%k]]);
                faces.push_back(q);
            }
            for (size_t i = 0; i != k; ++i)         // around the cell point
            {
                std::vector<Vertex> q;
                q.push_back(edgept[ce[i]]);
                q.push_back(facept[cf[(i+k-1)%k]]);
                q.push_back(C);
                q.push_back(facept[cf[i]]);
                faces.push_back(q);
            }

            const Cell nc2 = out.add_cell(faces);
            GISMO_ENSURE(nc2.is_valid(),
                         "gsVolCatmullClark: failed to build the cell at corner "
                         <<v<<" of cell "<<c<<".");
        }
    }

    *this->m_mesh = give(out);
}

} // namespace gismo
