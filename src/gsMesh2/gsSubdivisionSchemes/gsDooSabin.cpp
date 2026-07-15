/** @file gsDooSabin.cpp

    @brief Doo-Sabin subdivision on a  mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, D. Tolis, L. Mussmaecher
*/

#include <gsMesh2/gsSubdivisionSchemes/gsDooSabin.h>
#include <gsMesh2/gsSurfMesh.h>

namespace gismo
{

typedef gsSurfMesh::Point Point;
typedef gsSurfMesh::Vertex Vertex;
typedef gsSurfMesh::Face Face;
typedef gsSurfMesh::Halfedge Halfedge;
typedef gsSurfMesh::Edge Edge;

void gsDooSabin::subdivide_impl()
{
       
    index_t option = m_options.askInt("ds.boundaryMask");

    // New Mesh instance
    gsSurfMesh new_mesh;

    new_mesh.reserve(2 * m_mesh->n_edges(),
                     2 * (m_mesh->n_vertices() + m_mesh->n_faces() + m_mesh->n_edges()),
                     m_mesh->n_vertices() + m_mesh->n_faces() + m_mesh->n_edges());

    // Make a map to identify the new vertices
    std::map<Halfedge, Vertex> Map;
        
    // Create the vertices images and assign them to a halfedge
    for (auto he : m_mesh->halfedges()) // for all halfedges
    {
        if (!m_mesh->is_boundary(he))
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
    for (auto v : m_mesh->vertices()) 
    {
        ffv.clear();
        for (auto he : m_mesh->halfedges(v))
        {
            if (!m_mesh->is_boundary(he))
            {
                ffv.push_back(Map[he]);
            }
        }   

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
    for (auto f : m_mesh->faces())
    {
        ffv.clear();
        for (auto he : m_mesh->halfedges(f))
            ffv.push_back(Map[he]);
        new_mesh.add_face(ffv); // F-face
    }

    Halfedge h;
    // For all non-boundary edges create E-Faces
    for (auto ee : m_mesh->edges()) 
    {
        if (!(m_mesh->is_boundary(ee)))
        {
            ffv.clear();
            h = m_mesh->halfedge(ee, 0);
            ffv.push_back(Map[h]);
            h = m_mesh->cw_rotated_halfedge(h);
            ffv.push_back(Map[h]);
            h = m_mesh->prev_halfedge(h);
            ffv.push_back(Map[h]);
            h = m_mesh->cw_rotated_halfedge(h);
            ffv.push_back(Map[h]);
            new_mesh.add_face(ffv); // E-face
        }
    }

    *m_mesh = std::move(new_mesh);
}

gsSurfMesh::Point
gsDooSabin::ds_image_point_calc_interpolation(Halfedge hf)
{
    gsSurfMesh::Point coords;
    coords.setZero();
    Vertex oldv = m_mesh->from_vertex(hf);

    if (m_mesh->is_boundary(oldv)) 
    {

        if (m_mesh->valence(oldv) == 2) // Corner's case (same position because it will go to limit)
        {
            coords = m_mesh->position(oldv);
        }
        else // Chaikin method for the boundary
        { 
            if (m_mesh->is_boundary(m_mesh->to_vertex(hf))) 
            {
                coords = 0.75 * m_mesh->position(oldv) + 0.25 * m_mesh->position(m_mesh->to_vertex(hf));
            }
            else 
            {
                coords = 0.25 * m_mesh->position(m_mesh->from_vertex(m_mesh->prev_halfedge(hf))) + 0.75 * m_mesh->position(oldv);
            }
        }
    }
    else 
    {
        coords = ds_image_point_calc(hf);
    }
    return coords;
}


gsSurfMesh::Point
gsDooSabin::ds_image_point_calc(Halfedge hf)
{
    unsigned int face_valence{ m_mesh->valence(m_mesh->face(hf)) };

    Point coords;
    coords.setZero();

    int tempj = 0;
    real_t val = 0.0;
    // Coefficient of the first (current) vertex (i=j case)
    val = (real_t)(face_valence + 5) / (4 * face_valence);
    coords += val * m_mesh->position(m_mesh->from_vertex(hf));
    real_t sum_val = val;
    Halfedge next_he = m_mesh->next_halfedge(hf);

    // Creating the mask (image vertice coefficients) by looping to the number
    // of halfedges until I reach the initial half-edge (case i!=j).
    while (next_he != hf)
    {
        tempj++;
        val = (3 + 2 * cos(2.0 * EIGEN_PI * (0 - tempj) / face_valence)) / (4 * face_valence);
        coords += val * m_mesh->position(m_mesh->from_vertex(next_he));
        next_he = m_mesh->next_halfedge(next_he);
        sum_val += val;

    }

    return coords;

}

gsSurfMesh::Face_property<Point> gsDooSabin::face_limits(std::string label)
{

    Halfedge hb;
    auto limits = m_mesh->add_face_property<Point>((label), Point(0, 0, 0));
#   pragma omp parallel for default(shared)
    for (auto fit = m_mesh->faces_begin(); fit < m_mesh->faces_end(); ++fit)
    {
        if (m_mesh->is_boundary(*fit)) // Chaikin: middle edge point in the boundary
        {
            for (auto hit : m_mesh->halfedges(*fit))
                if (m_mesh->touches_boundary(hit))
                    hb = hit;
            if (m_mesh->touches_boundary(m_mesh->next_halfedge(hb))) // Corner case
                limits[*fit] = m_mesh->position(m_mesh->to_vertex(hb));
            else if (m_mesh->valence(m_mesh->from_vertex(hb)) == 2) // Corner case
                limits[*fit] = m_mesh->position(m_mesh->from_vertex(hb));
            else // Normal boundary
            {
                limits[*fit] = 0.5 * (m_mesh->position(m_mesh->from_vertex(hb)) +
                    m_mesh->position(m_mesh->to_vertex(hb)));
            }

            continue;
        }
        limits[*fit] = m_mesh->face_barycenter(*fit);
    }

    return limits;
}

gsSurfMesh::Face_property<Point> gsDooSabin::face_normal_limits(std::string label)
{
    bool normalize = m_options.askSwitch("normalize");
    auto points = m_mesh->points();
    auto limits = m_mesh->add_face_property<Point>((label), Point(0, 0, 0));
    unsigned int n;
    Halfedge he, hh;
    int i;
    Point t1, t2;
#   pragma omp parallel for default(shared) private(n,he,hh,i,t1,t2)
    for (auto fit = m_mesh->faces_begin(); fit < m_mesh->faces_end(); ++fit)
    {
        n = m_mesh->valence(*fit);
        if (m_mesh->is_boundary(*fit)) // TODO: Deprecate this feature by using middle edge point
        {
            gsWarn << "Boundary face is ignored.\n";
            continue;
        }

        auto& pt = limits[*fit];
        he = m_mesh->halfedge(*fit);
        hh = he;
        t1.setZero();
        t2.setZero();
        i = 0;
        do
        {
            t1 += math::cos(2 * i * EIGEN_PI / n) * points[m_mesh->from_vertex(hh)];
            t2 += math::sin(2 * i * EIGEN_PI / n) * points[m_mesh->from_vertex(hh)];
            hh = m_mesh->next_halfedge(hh);
            i++;
        } while (hh != he);
        pt = t1.cross(t2);

        if (normalize)
            pt = pt.normalized();
    }

    return limits;
}

} // namespace gismo
