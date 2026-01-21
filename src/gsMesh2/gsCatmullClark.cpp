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

    typedef gsSurfMesh::Point Point;
    typedef gsSurfMesh::Vertex Vertex;
    typedef gsSurfMesh::Face Face;
    typedef gsSurfMesh::Halfedge Halfedge;
    typedef gsSurfMesh::Edge Edge;

    void gsCatmullClark::subdivide(gsSurfMesh* mesh) {
        gsSurfMesh* m_mesh = mesh;
      
        gsSurfMesh::Vertex v;
        gsSurfMesh::Halfedge he;

        // reserve vertices, edges, faces
        m_mesh->reserve(m_mesh->n_vertices() + m_mesh->n_edges() + m_mesh->n_faces(),
            2 * m_mesh->n_edges(), 4 * m_mesh->n_faces());
        std::vector<std::string> he_props = m_mesh->halfedge_properties();
        const bool do_sharp = std::find(he_props.begin(), he_props.end(), "h:sharp") != he_props.end();
        gsSurfMesh::Halfedge_property<bool> sharp;
        if (do_sharp)
            sharp = m_mesh->get_halfedge_property<bool>("h:sharp");

        index_t env = m_mesh->n_vertices(); // edge vertices start here

        // loop over all edges, add edge points
        Point tmp;
        for (auto eit : m_mesh->edges())
        {
            he = m_mesh->halfedge(eit, 0);
            if (!(m_mesh->is_boundary(he) && m_mesh->is_boundary(m_mesh->opposite_halfedge(he))))
            {
                tmp = (m_mesh->position(m_mesh->from_vertex(he)) + m_mesh->position(m_mesh->to_vertex(he))) / 2;
                v = m_mesh->add_vertex(tmp);// edge points initialized as midpoints
                m_mesh->insert_vertex(he, v);
            }
        }

        index_t fnv = m_mesh->n_vertices(); // face vertices start here

        // loop over all faces, add face points
        for (auto fit : m_mesh->faces())
        {
            auto fv = m_mesh->vertices(fit);
            tmp.setZero();
            for (auto vc = fv.begin(); vc != fv.end(); ++vc, ++vc)
                tmp += m_mesh->position(*vc);
            tmp /= 4;
            m_mesh->add_vertex(tmp);  // vertex gets shifted face id
        }

        int i = 0;
        for (auto fit : m_mesh->faces())
        {
            v = gsSurfMesh::Vertex(fnv + (i++));//face vertex id ?
            //Start from an original vertex
            auto fv = m_mesh->vertices(fit).begin();
            if ((*fv).idx() >= env) ++fv; //todo: add -> operator
            //assert ( (*fv).idx() < nv )
            m_mesh->quad_split(fit, v, fv.he());
        }

#   pragma omp parallel for default(shared) private(v,tmp)
        for (i = env; i < fnv; ++i)
        {
            v = gsSurfMesh::Vertex(i); //edge points

            if (do_sharp && m_mesh->has_flag(v, sharp))
                continue; //remain midpoints
            else if (m_mesh->is_boundary(v))
            {
                //gsWarn<< "Boundary vertex "<< v.idx() <<"\n";
                continue;
            }
            else
            {
                //Interior edge points finalized here
                auto vit = m_mesh->vertices(v);
                auto vcp = vit;
                tmp.setZero();
                if (vit) do
                {
                    tmp += m_mesh->position(*vit);
                } while (++vit != vcp);
                tmp /= 4; // =valence(v);
                m_mesh->position(v) = tmp;
            }
        }

#   pragma omp parallel for default(shared) private(v)
        for (i = 0; i < env; ++i)
        {
            v = gsSurfMesh::Vertex(i); // original vertices
            auto n = m_mesh->valence(v);

            // Rigid corner ?
            if (n <= 2) continue;

            //original vertex positions are computed using new edge/face points only
            auto& pt = m_mesh->position(v);

            if (do_sharp && m_mesh->has_flag(v, sharp))
            {
                unsigned int sd = m_mesh->hcount(v, sharp);
                if (2 == sd)
                {
                    auto h1 = m_mesh->halfedge(v);
                    while (!sharp[h1] || m_mesh->touches_boundary(h1))
                        h1 = m_mesh->cw_rotated_halfedge(h1);
                    Halfedge h2 = h1;
                    h1 = m_mesh->ccw_rotated_halfedge(m_mesh->halfedge(v));
                    while (!sharp[h1] || m_mesh->is_boundary(h1))
                        h1 = m_mesh->ccw_rotated_halfedge(h1);

                    pt *= 2; //2~6
                    pt += m_mesh->position(m_mesh->to_vertex(h2)); // first flagged neighbor
                    pt += m_mesh->position(m_mesh->to_vertex(h1)); // second flagged neighbor
                    pt /= 4; // /4~8
                }

                if (sd != 1) // dart: use smooth rule, else  no relocation
                    continue;
            }

            if (m_mesh->is_boundary(v))
            {
                pt *= 2;
                auto hh = m_mesh->halfedge(v);
                while (!(m_mesh->touches_boundary(hh)))
                    hh = m_mesh->cw_rotated_halfedge(hh);
                //GISMO_ASSERT(touches_boundary(hh), "Did not find a boundary halfedge.");
                pt += m_mesh->position(m_mesh->to_vertex(hh)); //  right boundary neighbor
                hh = m_mesh->next_halfedge(m_mesh->opposite_halfedge(hh));
                pt += m_mesh->position(m_mesh->to_vertex(hh)); // left boundary neighbor
                pt /= 4;
            }
            else
            {
                auto vit = m_mesh->halfedges(v);
                auto vcp = vit;
                //formula: pt = ( (n*(n-3))*points[v] + 4*E - F ) / (n*n);
                pt *= n * (n - 3);
                if (vit)
                    do
                    { //pt += 4*E-F
                        pt += 4 * m_mesh->position(m_mesh->to_vertex(*vit))
                            - m_mesh->position(m_mesh->to_vertex(m_mesh->next_halfedge(*vit)));
                    } while (++vit != vcp);
                pt /= n * n;
            }
        }
    }

}
