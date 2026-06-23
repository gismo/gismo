/** @file gsCatmullClark.cpp

    @brief Catmull-Clark subdivision on a  mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, D. Tolis, L. Mussmaecher
*/


#include <gsMesh2/gsCatmullClark.h>
#include <gsMesh2/gsSubdivScheme.h>

namespace gismo {

void gsCatmullClark::subdivide_impl()
{
    GISMO_ASSERT(nullptr!=m_mesh,"Invalid mesh");
    gsCatmullClark::apply(*m_mesh);
}

void gsCatmullClark::apply(gsSurfMesh& mesh)
{
    gsSurfMesh::Vertex v;
    gsSurfMesh::Halfedge he;

    // reserve vertices, edges, faces
    mesh.reserve(mesh.n_vertices() + mesh.n_edges() + mesh.n_faces(),
                    2 * mesh.n_edges(), 4 * mesh.n_faces());
    std::vector<std::string> he_props = mesh.halfedge_properties();
    const bool do_sharp = std::find(he_props.begin(), he_props.end(), "h:sharp") != he_props.end();
    gsSurfMesh::Halfedge_property<bool> sharp;
    if (do_sharp)
        sharp = mesh.get_halfedge_property<bool>("h:sharp");

    index_t env = mesh.n_vertices(); // edge vertices start here

    // loop over all edges, add edge points
    Point tmp;
    for (auto eit : mesh.edges())
    {
        he = mesh.halfedge(eit, 0);
        if (!(mesh.is_boundary(he) && mesh.is_boundary(mesh.opposite_halfedge(he))))
        {
            tmp = (mesh.position(mesh.from_vertex(he)) + mesh.position(mesh.to_vertex(he))) / 2;
            v = mesh.add_vertex(tmp);// edge points initialized as midpoints
            mesh.insert_vertex(he, v);
        }
    }

    index_t fnv = mesh.n_vertices(); // face vertices start here

    // loop over all faces, add face points
    for (auto fit : mesh.faces())
    {
        auto fv = mesh.vertices(fit);
        tmp.setZero();
        for (auto vc = fv.begin(); vc != fv.end(); ++vc, ++vc)
            tmp += mesh.position(*vc);
        tmp /= 4;
        mesh.add_vertex(tmp);  // vertex gets shifted face id
    }

    int i = 0;
    for (auto fit : mesh.faces())
    {
        v = gsSurfMesh::Vertex(fnv + (i++));//face vertex id ?
        //Start from an original vertex
        auto fv = mesh.vertices(fit).begin();
        if ((*fv).idx() >= env) ++fv; //todo: add -> operator
        //assert ( (*fv).idx() < nv )
        mesh.quad_split(fit, v, fv.he());
    }

#   pragma omp parallel for default(shared) private(v,tmp)
    for (i = env; i < fnv; ++i)
    {
        v = gsSurfMesh::Vertex(i); //edge points

        if (do_sharp && mesh.has_flag(v, sharp))
            continue; //remain midpoints
        else if (mesh.is_boundary(v))
        {
            //gsWarn<< "Boundary vertex "<< v.idx() <<"\n";
            continue;
        }
        else
        {
            //Interior edge points finalized here
            auto vit = mesh.vertices(v);
            auto vcp = vit;
            tmp.setZero();
            if (vit) do
                     {
                         tmp += mesh.position(*vit);
                     } while (++vit != vcp);
            tmp /= 4; // =valence(v);
            mesh.position(v) = tmp;
        }
    }

#   pragma omp parallel for default(shared) private(v)
    for (i = 0; i < env; ++i)
    {
        v = gsSurfMesh::Vertex(i); // original vertices
        auto n = mesh.valence(v);

        // Rigid corner ?
        if (n <= 2) continue;

        //original vertex positions are computed using new edge/face points only
        auto& pt = mesh.position(v);

        if (do_sharp && mesh.has_flag(v, sharp))
        {
            unsigned int sd = mesh.hcount(v, sharp);
            if (2 == sd)
            {
                auto h1 = mesh.halfedge(v);
                while (!sharp[h1] || mesh.touches_boundary(h1))
                    h1 = mesh.cw_rotated_halfedge(h1);
                Halfedge h2 = h1;
                h1 = mesh.ccw_rotated_halfedge(mesh.halfedge(v));
                while (!sharp[h1] || mesh.is_boundary(h1))
                    h1 = mesh.ccw_rotated_halfedge(h1);

                pt *= 2; //2~6
                pt += mesh.position(mesh.to_vertex(h2)); // first flagged neighbor
                pt += mesh.position(mesh.to_vertex(h1)); // second flagged neighbor
                pt /= 4; // /4~8
            }

            if (sd != 1) // dart: use smooth rule, else  no relocation
                continue;
        }

        if (mesh.is_boundary(v))
        {
            pt *= 2;
            auto hh = mesh.halfedge(v);
            while (!(mesh.touches_boundary(hh)))
                hh = mesh.cw_rotated_halfedge(hh);
            //GISMO_ASSERT(touches_boundary(hh), "Did not find a boundary halfedge.");
            pt += mesh.position(mesh.to_vertex(hh)); //  right boundary neighbor
            hh = mesh.next_halfedge(mesh.opposite_halfedge(hh));
            pt += mesh.position(mesh.to_vertex(hh)); // left boundary neighbor
            pt /= 4;
        }
        else
        {
            auto vit = mesh.halfedges(v);
            auto vcp = vit;
            //formula: pt = ( (n*(n-3))*points[v] + 4*E - F ) / (n*n);
            pt *= n * (n - 3);
            if (vit)
                do
                { //pt += 4*E-F
                    pt += 4 * mesh.position(mesh.to_vertex(*vit))
                        - mesh.position(mesh.to_vertex(mesh.next_halfedge(*vit)));
                } while (++vit != vcp);
            pt /= n * n;
        }
    }
}

gsSurfMesh::Vertex_property<gsSurfMesh::Point> gsCatmullClark::vertex_limits(std::string label)
{
    auto points = m_mesh->points();
    auto limits = m_mesh->add_vertex_property<Point>(
        (label == "v:point" ? "v:limit_points_2022" : label), Point(0, 0, 0));
    real_t n;
#   pragma omp parallel for default(shared) private(n)
    for (auto vit = m_mesh->vertices_begin(); vit < m_mesh->vertices_end(); ++vit)
    {
        n = m_mesh->valence(*vit);
        if (m_mesh->is_boundary(*vit))
        {
            gsWarn << "Boundary vertex is ignored.\n";

            if (2 > n)
            {

            }
            continue;
        }

        auto& pt = limits[*vit];
        pt = n * n * points[*vit];
        for (auto he : m_mesh->halfedges(*vit))
        {
            if (m_mesh->is_boundary(he))
            {
                gsWarn << "Boundary halfedge is ignored.\n";
            }

            pt += 4 * points[m_mesh->to_vertex(he)] +
                points[m_mesh->to_vertex(m_mesh->next_halfedge(he))];
        }
        pt /= (n * (n + 5));
    }

    if (label == "v:point") //vertices are replaced by their limit positions
    {
        m_mesh->rename_vertex_property(points, "v:point_original");
        m_mesh->rename_vertex_property(limits, "v:point");
    }
    return limits;
}

gsSurfMesh::Vertex_property<gsSurfMesh::Point> gsCatmullClark::vertex_normal_limits(std::string label, bool normalize)
{
    auto points = m_mesh->points();
    //todo: check if label exists
    auto limits = m_mesh->add_vertex_property<Point>(label, Point(0, 0, 0));
    Point t1, t2;
    real_t c1, c2, cc1, cc2;
    index_t i;
    gsSurfMesh::Halfedge h2;
#   pragma omp parallel for default(shared) private(h2,t1,t2,c1,c2,cc1,cc2,i)
    for (auto vit = m_mesh->vertices_begin(); vit < m_mesh->vertices_end(); ++vit)
    {
        const real_t n = m_mesh->valence(*vit);
        const real_t cospin = math::cos(EIGEN_PI / n);
        cc2 = 1 / (n * math::sqrt(4 + cospin * cospin));
        cc1 = 1 / n + cospin * cc2;
        t1.setZero();
        t2.setZero();
        i = 0;
        for (auto he : m_mesh->halfedges(*vit))
        {
            h2 = m_mesh->ccw_rotated_halfedge(he);
            c1 = math::cos(2 * i * EIGEN_PI / n) * cc1;
            c2 = math::cos((2 * i + 1) * EIGEN_PI / n) * cc2;
            t1 += c1 * points[m_mesh->to_vertex(he)]
                + c2 * points[m_mesh->to_vertex(m_mesh->next_halfedge(he))];
            t2 += c1 * points[m_mesh->to_vertex(h2)]
                + c2 * points[m_mesh->to_vertex(m_mesh->next_halfedge(h2))];
            ++i;
        }
        if (normalize)
            limits[*vit] = t1.cross(t2).normalized();
        else
            limits[*vit] = t1.cross(t2);
    }
    return limits;
}

gsSurfMesh::Vertex_property<gsSurfMesh::Point> gsCatmullClark::vertex_tangent_limits(std::string label, bool normalize)
{
    gsSurfMesh::Vertex v;
    gsSurfMesh::Halfedge h2;

    auto points = m_mesh->points();
    //todo: check if label exists
    auto limits = m_mesh->add_vertex_property<Point>(label, Point(0, 0, 0));
    Point t1, t2;
    real_t c1, c2, cc1, cc2;
    index_t i;
#   pragma omp parallel for default(shared) private(v,h2,t1,t2,c1,c2,cc1,cc2,i)
    for (auto vit = m_mesh->vertices_begin(); vit < m_mesh->vertices_end(); ++vit)
    {
        const real_t n = m_mesh->valence(*vit);
        const real_t cospin = math::cos(EIGEN_PI / n);
        cc2 = 1 / (n * math::sqrt(4 + cospin * cospin));
        cc1 = 1 / n + cospin * cc2;
        t1.setZero();
        t2.setZero();
        i = 0;
        for (auto he : m_mesh->halfedges(*vit))
        {
            h2 = m_mesh->ccw_rotated_halfedge(he);
            c1 = math::cos(2 * i * EIGEN_PI / n) * cc1;
            c2 = math::cos((2 * i + 1) * EIGEN_PI / n) * cc2;
            t1 += c1 * points[m_mesh->to_vertex(he)]
                + c2 * points[m_mesh->to_vertex(m_mesh->next_halfedge(he))];
            ++i;
        }
        if (normalize)
            limits[*vit] = t1.normalized();
        else
            limits[*vit] = t1;
    }
    return limits;
}


}
