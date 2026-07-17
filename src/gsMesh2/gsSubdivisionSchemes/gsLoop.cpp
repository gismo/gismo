/** @file gsLoop.cpp

    @brief Loop subdivision on a triangular mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, D. Tolis, L. Mussmaecher
*/

#include <gsMesh2/gsSubdivisionSchemes/gsLoop.h>
#include <gsMesh2/gsSurfMesh.h>

namespace gismo
{

typedef gsSurfMesh::Point Point;
typedef gsSurfMesh::Vertex Vertex;
typedef gsSurfMesh::Face Face;
typedef gsSurfMesh::Halfedge Halfedge;
typedef gsSurfMesh::Edge Edge;

void gsLoop::subdivide_impl()
{
    index_t mask_option = m_options.askInt("loop.maskType");

    gsSurfMesh nm;
    Vertex v;

    // Tessalate the current mesh, in order to have pure triangle mesh
    for (auto fit : m_mesh->faces())
        if (m_mesh->valence(fit) > 3)
            m_mesh->triangulate(fit);


    // reserve vertices, edges, faces
    m_mesh->reserve(m_mesh->n_vertices() + m_mesh->n_edges() + m_mesh->n_faces(),
        2 * m_mesh->n_edges(), 4 * m_mesh->n_faces());


    index_t fnv = m_mesh->n_vertices(); // Number of vertices in the current mesh
    index_t num_f = m_mesh->n_faces(); // Number of faces in the current mesh

    gsSurfMesh::Halfedge he, nh, ph, hh;
    Point tmp;
    std::map<Edge, Vertex> edgeVerts;
    std::map<Vertex, Vertex> innerVerts;
    // In each edge create the new vertices  (odd rules)
    for (auto ho : m_mesh->halfedges())
    {

        if (m_mesh->touches_boundary(ho))
        {
            tmp = (real_t)1 / 2 * m_mesh->position(m_mesh->from_vertex(ho)) +
                (real_t)1 / 2 * m_mesh->position(m_mesh->to_vertex(ho));
        }
        else if (m_mesh->is_boundary(ho))
        {
            continue;
        }
        else
        {
            nh = m_mesh->next_halfedge(ho);
            ph = m_mesh->prev_halfedge(m_mesh->opposite_halfedge(ho));
            tmp = (real_t)3 / 8 * m_mesh->position(m_mesh->from_vertex(ho)) +
                (real_t)3 / 8 * m_mesh->position(m_mesh->to_vertex(ho)) +
                (real_t)1 / 8 * m_mesh->position(m_mesh->to_vertex(nh)) +
                (real_t)1 / 8 * m_mesh->position(m_mesh->from_vertex(ph));
        }
        if (edgeVerts.count(m_mesh->edge(ho)) > 0)
        {
            continue;
        }
        v = nm.add_vertex(tmp);

        edgeVerts[m_mesh->edge(ho)] = v;
    }

    index_t n;
    real_t a, b, c;

    // For each old vertex modify its position (even rules)
    // For boundary we will impose the classic rule which causes cubic spline at the infinity
    // In case of boundary vertices in corner we will leave them in tact (See Hoppe 1994, Bierman 2000, Ling 2008)
    for (auto it = 0; it < fnv; it++)
    {
        v = Vertex(it);
        if (m_mesh->is_boundary(v))
        {


            he = m_mesh->halfedge(v);
            hh = he;

            // Find he touching boundary
            while (!(m_mesh->touches_boundary(hh)))
            {
                hh = m_mesh->ccw_rotated_halfedge(hh);
            }
            he = hh;
            // Apply boundary mask
            tmp = (real_t)3 / 4 * m_mesh->position(m_mesh->from_vertex(he)) +
                (real_t)1 / 8 * m_mesh->position(m_mesh->to_vertex(he));
            hh = he;
            do
            {
                hh = m_mesh->ccw_rotated_halfedge(hh);
                if (m_mesh->is_boundary(m_mesh->to_vertex(hh)))
                {
                    tmp += (real_t)1 / 8 * m_mesh->position(m_mesh->to_vertex(hh));
                }
            } while (m_mesh->ccw_rotated_halfedge(hh) != he);
        }
        else
        {
            n = m_mesh->valence(v);
            he = m_mesh->halfedge(v);
            if (mask_option == 0)
            {  // Simplified averganging mask (Warren,Weimer 2002)
                if (n == 3)
                {
                    a = (real_t)3 / 16;
                }
                else
                {
                    a = (real_t)3 / (8 * n);
                }
            }
            else
            {  // Loop's original mask (1987)
                a = (0.625 - (0.375 + 0.25 * cos(2 * 3.14159 / (real_t)n)) * (0.375 + 0.25 * cos(2 * 3.14159 / (real_t)n))) / (real_t)n;
            }
            b = a;
            c = 1.0 - a * n;



            // Apply inner mask
            tmp = c * m_mesh->position(v);
            hh = he;
            do
            {
                tmp += b * m_mesh->position(m_mesh->to_vertex(hh));
                hh = m_mesh->ccw_rotated_halfedge(hh);
            } while (hh != he);

        }
        innerVerts[v] = nm.add_vertex(tmp);
    }

    for (auto h : m_mesh->halfedges())
    {

        if (m_mesh->is_boundary(h))
        {
            continue;
        }

        v = m_mesh->from_vertex(h);
        ph = m_mesh->prev_halfedge(h);
        nm.add_triangle(innerVerts[v], edgeVerts[m_mesh->edge(h)], edgeVerts[m_mesh->edge(ph)]);
    }
    Edge e;
    Face f;
    std::vector<Vertex> edgeVertsFace;
    for (auto it = 0; it < num_f; it++)
    {
        f = Face(it);
        edgeVertsFace.clear();
        for (auto hit : m_mesh->halfedges(f))
        {
            e = m_mesh->edge(hit);
            edgeVertsFace.push_back(edgeVerts[e]);
        }

        nm.add_face(edgeVertsFace);

    }


    *m_mesh = std::move(nm);
}


} // namespace gismo
