/** @file gsLoop.hpp

    @brief Loop subdivision on a triangular mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, D. Tolis, L. Mussmaecher
*/

#pragma once

#include <gsMesh2/gsSubdivisionSchemes/gsLoop.h>
#include <gsMesh2/gsSurfMesh.h>

namespace gismo
{

template <class Scalar>
void gsLoop<Scalar>::subdivide_impl()
{
    typedef typename gsSurfMesh<Scalar>::Point Point;
    typedef typename gsSurfMesh<Scalar>::Vertex Vertex;
    typedef typename gsSurfMesh<Scalar>::Face Face;
    typedef typename gsSurfMesh<Scalar>::Edge Edge;
    typedef typename gsSurfMesh<Scalar>::Halfedge Halfedge;

    index_t mask_option = this->m_options.askInt("loop.maskType");

    gsSurfMesh<Scalar> nm;
    Vertex v;

    // Tessalate the current mesh, in order to have pure triangle mesh
    for (auto fit : this->m_mesh->faces())
        if (this->m_mesh->valence(fit) > 3)
            this->m_mesh->triangulate(fit);


    // reserve vertices, edges, faces
    this->m_mesh->reserve(this->m_mesh->n_vertices() + this->m_mesh->n_edges() + this->m_mesh->n_faces(),
        2 * this->m_mesh->n_edges(), 4 * this->m_mesh->n_faces());


    index_t fnv = this->m_mesh->n_vertices(); // Number of vertices in the current mesh
    index_t num_f = this->m_mesh->n_faces(); // Number of faces in the current mesh

    Halfedge he, nh, ph, hh;
    Point tmp;
    std::map<Edge, Vertex> edgeVerts;
    std::map<Vertex, Vertex> innerVerts;
    // In each edge create the new vertices  (odd rules)
    for (auto ho : this->m_mesh->halfedges())
    {

        if (this->m_mesh->touches_boundary(ho))
        {
            tmp = (Scalar)1 / 2 * this->m_mesh->position(this->m_mesh->from_vertex(ho)) +
                (Scalar)1 / 2 * this->m_mesh->position(this->m_mesh->to_vertex(ho));
        }
        else if (this->m_mesh->is_boundary(ho))
        {
            continue;
        }
        else
        {
            nh = this->m_mesh->next_halfedge(ho);
            ph = this->m_mesh->prev_halfedge(this->m_mesh->opposite_halfedge(ho));
            tmp = (Scalar)3 / 8 * this->m_mesh->position(this->m_mesh->from_vertex(ho)) +
                (Scalar)3 / 8 * this->m_mesh->position(this->m_mesh->to_vertex(ho)) +
                (Scalar)1 / 8 * this->m_mesh->position(this->m_mesh->to_vertex(nh)) +
                (Scalar)1 / 8 * this->m_mesh->position(this->m_mesh->from_vertex(ph));
        }
        if (edgeVerts.count(this->m_mesh->edge(ho)) > 0)
        {
            continue;
        }
        v = nm.add_vertex(tmp);

        edgeVerts[this->m_mesh->edge(ho)] = v;
    }

    index_t n;
    Scalar a, b, c;

    // For each old vertex modify its position (even rules)
    // For boundary we will impose the classic rule which causes cubic spline at the infinity
    // In case of boundary vertices in corner we will leave them in tact (See Hoppe 1994, Bierman 2000, Ling 2008)
    for (auto it = 0; it < fnv; it++)
    {
        v = Vertex(it);
        if (this->m_mesh->is_boundary(v))
        {


            he = this->m_mesh->halfedge(v);
            hh = he;

            // Find he touching boundary
            while (!(this->m_mesh->touches_boundary(hh)))
            {
                hh = this->m_mesh->ccw_rotated_halfedge(hh);
            }
            he = hh;
            // Apply boundary mask
            tmp = (Scalar)3 / 4 * this->m_mesh->position(this->m_mesh->from_vertex(he)) +
                (Scalar)1 / 8 * this->m_mesh->position(this->m_mesh->to_vertex(he));
            hh = he;
            do
            {
                hh = this->m_mesh->ccw_rotated_halfedge(hh);
                if (this->m_mesh->is_boundary(this->m_mesh->to_vertex(hh)))
                {
                    tmp += (Scalar)1 / 8 * this->m_mesh->position(this->m_mesh->to_vertex(hh));
                }
            } while (this->m_mesh->ccw_rotated_halfedge(hh) != he);
        }
        else
        {
            n = this->m_mesh->valence(v);
            he = this->m_mesh->halfedge(v);
            if (mask_option == 0)
            {  // Simplified averganging mask (Warren,Weimer 2002)
                if (n == 3)
                {
                    a = (Scalar)3 / 16;
                }
                else
                {
                    a = (Scalar)3 / (8 * n);
                }
            }
            else
            {  // Loop's original mask (1987)
                a = (0.625 - (0.375 + 0.25 * math::cos(2 * 3.14159 / (Scalar)n)) * (0.375 + 0.25 * math::cos(2 * 3.14159 / (Scalar)n))) / (Scalar)n;
            }
            b = a;
            c = 1.0 - a * n;



            // Apply inner mask
            tmp = c * this->m_mesh->position(v);
            hh = he;
            do
            {
                tmp += b * this->m_mesh->position(this->m_mesh->to_vertex(hh));
                hh = this->m_mesh->ccw_rotated_halfedge(hh);
            } while (hh != he);

        }
        innerVerts[v] = nm.add_vertex(tmp);
    }

    for (auto h : this->m_mesh->halfedges())
    {

        if (this->m_mesh->is_boundary(h))
        {
            continue;
        }

        v = this->m_mesh->from_vertex(h);
        ph = this->m_mesh->prev_halfedge(h);
        nm.add_triangle(innerVerts[v], edgeVerts[this->m_mesh->edge(h)], edgeVerts[this->m_mesh->edge(ph)]);
    }
    Edge e;
    Face f;
    std::vector<Vertex> edgeVertsFace;
    for (auto it = 0; it < num_f; it++)
    {
        f = Face(it);
        edgeVertsFace.clear();
        for (auto hit : this->m_mesh->halfedges(f))
        {
            e = this->m_mesh->edge(hit);
            edgeVertsFace.push_back(edgeVerts[e]);
        }

        nm.add_face(edgeVertsFace);

    }


    *this->m_mesh = std::move(nm);
}


} // namespace gismo
