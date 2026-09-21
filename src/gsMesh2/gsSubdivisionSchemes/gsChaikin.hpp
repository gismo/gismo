/** @file gsChaikin.hpp

    @brief Chaikin subdivision on a polyline.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, D. Tolis, L. Mussmaecher
*/

#pragma once

#include <gsMesh2/gsSubdivisionSchemes/gsChaikin.h>
#include <gsMesh2/gsSurfMesh.h>

namespace gismo
{

template <class Scalar>
void gsChaikin<Scalar>::subdivide_impl()
{
    typedef typename gsSurfMesh<Scalar>::Point Point;
    typedef typename gsSurfMesh<Scalar>::Vertex Vertex;
    typedef typename gsSurfMesh<Scalar>::Halfedge Halfedge;

    // Splitting phase (adding new vertices at edges midpoints)
    Point tmp;
    Halfedge he;
    for (auto eit : this->m_mesh->edges())
    {
        he = this->m_mesh->halfedge(eit, 0);
        tmp = 0.5 * (this->m_mesh->position(this->m_mesh->from_vertex(he)) + this->m_mesh->position(this->m_mesh->to_vertex(he)));
        this->m_mesh->insert_vertex(eit, tmp);
    }

    he = this->m_mesh->halfedge(typename gsSurfMesh<Scalar>::Face(0)); // Face's halfedge (we consider closed curves for now)

    Halfedge hh = he;
    gsVector<Point> newverts;
    newverts.resize(this->m_mesh->n_vertices());
    index_t i = 0;

    // Averanging phase (Chaikin's scheme does the averenging in 0.5 of two neighbour vertices)
    do
    {
        tmp = 0.5 * (this->m_mesh->position(this->m_mesh->from_vertex(hh)) + this->m_mesh->position(this->m_mesh->to_vertex(hh)));
        newverts(i++) = tmp;
        hh = this->m_mesh->next_halfedge(hh);
    } while (he != hh);

    i = 0;
    hh = he;

    // Generating new curve
    do
    {
        this->m_mesh->position(this->m_mesh->to_vertex(hh)) = newverts(i++);
        hh = this->m_mesh->next_halfedge(hh);
    } while (he != hh);
}


} // namespace gismo
