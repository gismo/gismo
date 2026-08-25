/** @file gsDooSabin.hpp

    @brief Doo-Sabin subdivision on a  mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, D. Tolis, L. Mussmaecher
*/

#pragma once

#include <gsMesh2/gsSubdivisionSchemes/gsDooSabin.h>
#include <gsMesh2/gsSurfMesh.h>

namespace gismo
{

template <class Scalar>
void gsDooSabin<Scalar>::subdivide_impl()
{
    typedef typename gsSurfMesh<Scalar>::Vertex Vertex;
    typedef typename gsSurfMesh<Scalar>::Face Face;
    typedef typename gsSurfMesh<Scalar>::Halfedge Halfedge;

    index_t option = this->m_options.askInt("ds.boundaryMask");

    // New Mesh instance
    gsSurfMesh<Scalar> new_mesh;

    new_mesh.reserve(2 * this->m_mesh->n_edges(),
                     2 * (this->m_mesh->n_vertices() + this->m_mesh->n_faces() + this->m_mesh->n_edges()),
                     this->m_mesh->n_vertices() + this->m_mesh->n_faces() + this->m_mesh->n_edges());

    // Make a map to identify the new vertices
    std::map<Halfedge, Vertex> Map;

    // Create the vertices images and assign them to a halfedge
    for (auto he : this->m_mesh->halfedges()) // for all halfedges
    {
        if (!this->m_mesh->is_boundary(he))
        {
            if (option > 0) // trimmed cases (option = 1).
            {
                Map[he] = new_mesh.add_vertex(ds_image_point_calc(he));
            }
            else // interpolatory case (option = 0).
            {
                Map[he] = new_mesh.add_vertex(ds_image_point_calc_interpolation(he));
            }
        }
    }

    std::vector<Vertex> ffv;
    // For all vertices create V-faces
    for (auto v : this->m_mesh->vertices())
    {
        ffv.clear();
        for (auto he : this->m_mesh->halfedges(v))
        {
            if (!this->m_mesh->is_boundary(he))
            {
                ffv.push_back(Map[he]);
            }
        }

        if (ffv.size() < 2) // corner case
            continue;

        if (ffv.size() == 2)
        {
            new_mesh.add_edge(ffv[0], ffv[1]); // V-Edge in the regular boundary case
        }
        else
        {
            new_mesh.add_face(ffv); // V-face
        }
    }

    // For all faces create F-Faces
    for (auto f : this->m_mesh->faces())
    {
        ffv.clear();
        for (auto he : this->m_mesh->halfedges(f))
            ffv.push_back(Map[he]);
        new_mesh.add_face(ffv); // F-face
    }

    Halfedge h;
    // For all non-boundary edges create E-Faces
    for (auto ee : this->m_mesh->edges())
    {
        if (!(this->m_mesh->is_boundary(ee)))
        {
            ffv.clear();
            h = this->m_mesh->halfedge(ee, 0);
            ffv.push_back(Map[h]);
            h = this->m_mesh->cw_rotated_halfedge(h);
            ffv.push_back(Map[h]);
            h = this->m_mesh->prev_halfedge(h);
            ffv.push_back(Map[h]);
            h = this->m_mesh->cw_rotated_halfedge(h);
            ffv.push_back(Map[h]);
            new_mesh.add_face(ffv); // E-face
        }
    }

    *this->m_mesh = std::move(new_mesh);
}

template <class Scalar>
typename gsSurfMesh<Scalar>::Point
gsDooSabin<Scalar>::ds_image_point_calc_interpolation(Halfedge hf)
{
    typedef typename gsSurfMesh<Scalar>::Vertex Vertex;

    Point coords;
    coords.setZero();
    Vertex oldv = this->m_mesh->from_vertex(hf);

    if (this->m_mesh->is_boundary(oldv))
    {

        if (this->m_mesh->valence(oldv) == 2) // Corner's case (same position because it will go to limit)
        {
            coords = this->m_mesh->position(oldv);
        }
        else // Chaikin method for the boundary
        {
            if (this->m_mesh->is_boundary(this->m_mesh->to_vertex(hf)))
            {
                coords = 0.75 * this->m_mesh->position(oldv) + 0.25 * this->m_mesh->position(this->m_mesh->to_vertex(hf));
            }
            else
            {
                coords = 0.25 * this->m_mesh->position(this->m_mesh->from_vertex(this->m_mesh->prev_halfedge(hf))) + 0.75 * this->m_mesh->position(oldv);
            }
        }
    }
    else
    {
        coords = ds_image_point_calc(hf);
    }
    return coords;
}


template <class Scalar>
typename gsSurfMesh<Scalar>::Point
gsDooSabin<Scalar>::ds_image_point_calc(Halfedge hf)
{
    unsigned int face_valence{ this->m_mesh->valence(this->m_mesh->face(hf)) };

    Point coords;
    coords.setZero();

    int tempj = 0;
    Scalar val = 0.0;
    // Coefficient of the first (current) vertex (i=j case)
    val = (Scalar)(face_valence + 5) / (4 * face_valence);
    coords += val * this->m_mesh->position(this->m_mesh->from_vertex(hf));
    Scalar sum_val = val;
    Halfedge next_he = this->m_mesh->next_halfedge(hf);

    // Creating the mask (image vertice coefficients) by looping to the number
    // of halfedges until I reach the initial half-edge (case i!=j).
    while (next_he != hf)
    {
        tempj++;
        val = (3 + 2 * math::cos(2.0 * EIGEN_PI * (0 - tempj) / face_valence)) / (4 * face_valence);
        coords += val * this->m_mesh->position(this->m_mesh->from_vertex(next_he));
        next_he = this->m_mesh->next_halfedge(next_he);
        sum_val += val;

    }

    return coords;

}

template <class Scalar>
typename gsSurfMesh<Scalar>::template Face_property<typename gsSurfMesh<Scalar>::Point>
gsDooSabin<Scalar>::face_limits(std::string label)
{
    Halfedge hb;
    auto limits = this->m_mesh->template add_face_property<Point>((label), Point(0, 0, 0));
#   pragma omp parallel for default(shared)
    for (auto fit = this->m_mesh->faces_begin(); fit < this->m_mesh->faces_end(); ++fit)
    {
        if (this->m_mesh->is_boundary(*fit)) // Chaikin: middle edge point in the boundary
        {
            for (auto hit : this->m_mesh->halfedges(*fit))
                if (this->m_mesh->touches_boundary(hit))
                    hb = hit;
            if (this->m_mesh->touches_boundary(this->m_mesh->next_halfedge(hb))) // Corner case
                limits[*fit] = this->m_mesh->position(this->m_mesh->to_vertex(hb));
            else if (this->m_mesh->valence(this->m_mesh->from_vertex(hb)) == 2) // Corner case
                limits[*fit] = this->m_mesh->position(this->m_mesh->from_vertex(hb));
            else // Normal boundary
            {
                limits[*fit] = 0.5 * (this->m_mesh->position(this->m_mesh->from_vertex(hb)) +
                    this->m_mesh->position(this->m_mesh->to_vertex(hb)));
            }

            continue;
        }
        limits[*fit] = this->m_mesh->face_barycenter(*fit);
    }

    return limits;
}

template <class Scalar>
typename gsSurfMesh<Scalar>::template Face_property<typename gsSurfMesh<Scalar>::Point>
gsDooSabin<Scalar>::face_normal_limits(std::string label)
{
    bool normalize = this->m_options.askSwitch("normalize");
    auto points = this->m_mesh->points();
    auto limits = this->m_mesh->template add_face_property<Point>((label), Point(0, 0, 0));
    unsigned int n;
    Halfedge he, hh;
    int i;
    Point t1, t2;
#   pragma omp parallel for default(shared) private(n,he,hh,i,t1,t2)
    for (auto fit = this->m_mesh->faces_begin(); fit < this->m_mesh->faces_end(); ++fit)
    {
        n = this->m_mesh->valence(*fit);
        if (this->m_mesh->is_boundary(*fit)) // TODO: Deprecate this feature by using middle edge point
        {
            gsWarn << "Boundary face is ignored.\n";
            continue;
        }

        auto& pt = limits[*fit];
        he = this->m_mesh->halfedge(*fit);
        hh = he;
        t1.setZero();
        t2.setZero();
        i = 0;
        do
        {
            t1 += math::cos(2 * i * EIGEN_PI / n) * points[this->m_mesh->from_vertex(hh)];
            t2 += math::sin(2 * i * EIGEN_PI / n) * points[this->m_mesh->from_vertex(hh)];
            hh = this->m_mesh->next_halfedge(hh);
            i++;
        } while (hh != he);
        pt = t1.cross(t2);

        if (normalize)
            pt = pt.normalized();
    }

    return limits;
}

} // namespace gismo
