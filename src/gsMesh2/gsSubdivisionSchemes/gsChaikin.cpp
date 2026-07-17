/** @file gsDooSabin.cpp

    @brief Doo-Sabin subdivision on a  mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, D. Tolis, L. Mussmaecher
*/

#include <gsMesh2/gsSubdivisionSchemes/gsChaikin.h>
#include <gsMesh2/gsSurfMesh.h>

namespace gismo
{

typedef gsSurfMesh::Point Point;
typedef gsSurfMesh::Vertex Vertex;
typedef gsSurfMesh::Face Face;
typedef gsSurfMesh::Halfedge Halfedge;
typedef gsSurfMesh::Edge Edge;

void gsChaikin::subdivide_impl()
{
       
    // Splitting phase (adding new vertices at edges midpoints)
    gsSurfMesh::Point tmp, tmps, tmpe;
    gsSurfMesh::Halfedge he;
    for (auto eit : m_mesh->edges())
    {
        he = m_mesh->halfedge(eit, 0);
        tmp = 0.5 * (m_mesh->position(m_mesh->from_vertex(he)) + m_mesh->position(m_mesh->to_vertex(he)));
        m_mesh->insert_vertex(eit, tmp);
    }

    he = m_mesh->halfedge(gsSurfMesh::Face(0)); // Face's halfedge (we consider closed curves for now)

    gsSurfMesh::Halfedge hh = he;
    gsVector<gsSurfMesh::Point> newverts;
    newverts.resize(m_mesh->n_vertices());
    index_t i = 0;

    // Averanging phase (Chaikin's scheme does the averenging in 0.5 of two neighbour vertices)
    do
    {
        tmp = 0.5 * (m_mesh->position(m_mesh->from_vertex(hh)) + m_mesh->position(m_mesh->to_vertex(hh)));
        newverts(i++) = tmp;
        hh = m_mesh->next_halfedge(hh);
    } while (he != hh);

    i = 0;
    hh = he;

    // Generating new curve
    do
    {
        m_mesh->position(m_mesh->to_vertex(hh)) = newverts(i++);
        hh = m_mesh->next_halfedge(hh);
    } while (he != hh);
}


} // namespace gismo
