/** @file gsDooSabin.cpp

    @brief Doo-Sabin subdivision on a  mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, D. Tolis, L. Mussmaecher
*/


#include <gsMesh2/gsSurfMesh.h>
#include <gsMesh2/gsDooSabin.h>

namespace gismo {

    typedef gsSurfMesh::Point Point;
    typedef gsSurfMesh::Vertex Vertex;
    typedef gsSurfMesh::Face Face;
    typedef gsSurfMesh::Halfedge Halfedge;
    typedef gsSurfMesh::Edge Edge;

    void gsDooSabin::subdivide(gsSurfMesh* mesh) {
        gsSurfMesh* m_mesh = mesh;   
       
        index_t option = m_options.askInt("ds.boundaryMask");

        // New Mesh instance
        gsSurfMesh new_mesh;

        // Make a map to identify the new vertices
        std::map<std::pair<Vertex, Face>, Vertex> Map;


        // Create of V-Faces
        Vertex v;
        std::vector<Vertex> ffv;
        for (auto oldv : m_mesh->vertices()) {

            // Map vertices and face of old mesh with new (calculated) vertex
            ffv.clear();
            for (auto oldf : m_mesh->faces(oldv))
            {
                if (option == 1)
                    v = new_mesh.add_vertex(ds_image_point_calc_vanila(oldv, oldf, *m_mesh));
                else
                    v = new_mesh.add_vertex(ds_image_point_calc_interpolation(oldv, oldf, *m_mesh));

                ffv.push_back(v);
                Map[std::make_pair(oldv, oldf)] = v;
            }
            if (ffv.size() == 2) {
                new_mesh.add_edge(ffv[0], ffv[1]);

            }
            else if (ffv.size() > 2) {
                new_mesh.add_face(ffv);

            }
        }

        // Create of F-faces
        for (auto oldf : m_mesh->faces()) {
            ffv.clear();
            for (auto oldv : m_mesh->vertices(oldf)) {
                ffv.push_back(Map[std::make_pair(oldv, oldf)]);
            }
            new_mesh.add_face(ffv);
        }

        // Create E-Face by looping through boundary of F,V-faces
        for (auto olde : m_mesh->edges()) {
            if (!(m_mesh->is_boundary(olde))) {
                auto h0 = m_mesh->halfedge(olde, 0);
                auto h1 = m_mesh->halfedge(olde, 1);


                new_mesh.add_quad(Map[std::make_pair(m_mesh->from_vertex(h1), m_mesh->face(h0))],
                    Map[std::make_pair(m_mesh->from_vertex(h0), m_mesh->face(h0))],
                    Map[std::make_pair(m_mesh->from_vertex(h0), m_mesh->face(h1))],
                    Map[std::make_pair(m_mesh->from_vertex(h1), m_mesh->face(h1))]

                );
            }
        }
        *m_mesh = std::move(new_mesh);

    }

    gsSurfMesh::Point
    gsDooSabin::ds_image_point_calc_interpolation(Vertex oldv, Face oldf, const gsSurfMesh& mesh)
    {
        unsigned int face_valence{ mesh.valence(oldf) };

        // Find the halfedge of the vertex I am looking in case the IDs of mesh
        // are not sequencial (i.e. Quad ID: 1,23,3,5)
        Halfedge hf = mesh.halfedge(oldf);
        for (auto hh : mesh.halfedges(oldf))
        {
            if (mesh.from_vertex(hh) == oldv)
            {
                hf = hh;
                break;
            }
        }

        real_t val{ 0 };

        gsEigen::Matrix<double, 3, 1, 0, 3, 1> coords;
        coords.setZero();


        if (mesh.is_boundary(oldv)) {

            if (mesh.valence(oldv) == 2) { // Corner's case (same position because it will go to limit)
                coords = mesh.position(oldv);
            }
            else { // Chaikin method for the boundary
                if (mesh.is_boundary(mesh.to_vertex(hf))) {
                    coords = 0.75 * mesh.position(oldv) + 0.25 * mesh.position(mesh.to_vertex(hf));
                }
                else {
                    coords = 0.25 * mesh.position(mesh.from_vertex(mesh.prev_halfedge(hf))) + 0.75 * mesh.position(oldv);
                }
            }


        }
        else {

            int tempj{ 0 };

            // Coefficient of the first (current) vertex (i=j case)
            val = (real_t)(face_valence + 5) / (4 * face_valence);
            coords += val * mesh.position(mesh.from_vertex(hf));
            real_t sum_val{ val };
            Halfedge next_he{ mesh.next_halfedge(hf) };

            // Creating the mask (image vertice coefficients) by looping to the number
            // of halfedges until I reach the initial half-edge (case i!=j).
            while (next_he != hf) {
                tempj++;
                val = (3 + 2 * cos(2.0 * EIGEN_PI * (0 - tempj) / face_valence)) / (4 * face_valence);
                coords += val * mesh.position(mesh.from_vertex(next_he));
                next_he = mesh.next_halfedge(next_he);
                sum_val += val;

            }


        }


        Point temp;
        temp[0] = coords(0);
        temp[1] = coords(1);
        temp[2] = coords(2);

        return temp;

    }

    gsSurfMesh::Point
    gsDooSabin::ds_image_point_calc_vanila(Vertex oldv, Face oldf, const gsSurfMesh& mesh)
    {
        unsigned int face_valence{ mesh.valence(oldf) };

        // Find the halfedge of the vertex I am looking in case the IDs of mesh
        // are not sequencial (i.e. Quad ID: 1,23,3,5)
        Halfedge hf = mesh.halfedge(oldf);
        for (auto hh : mesh.halfedges(oldf))
        {
            if (mesh.from_vertex(hh) == oldv)
            {
                hf = hh;
                break;
            }
        }

        real_t val{ 0 };

        gsEigen::Matrix<double, 3, 1, 0, 3, 1> coords;
        coords.setZero();

        int tempj{ 0 };

        // Coefficient of the first (current) vertex (i=j case)
        val = (real_t)(face_valence + 5) / (4 * face_valence);
        coords += val * mesh.position(mesh.from_vertex(hf));
        real_t sum_val{ val };
        Halfedge next_he{ mesh.next_halfedge(hf) };

        // Creating the mask (image vertice coefficients) by looping to the number
        // of halfedges until I reach the initial half-edge (case i!=j).
        while (next_he != hf) {
            tempj++;
            val = (3 + 2 * cos(2.0 * EIGEN_PI * (0 - tempj) / face_valence)) / (4 * face_valence);
            coords += val * mesh.position(mesh.from_vertex(next_he));
            next_he = mesh.next_halfedge(next_he);
            sum_val += val;

        }

        Point temp;
        temp[0] = coords(0);
        temp[1] = coords(1);
        temp[2] = coords(2);

        return temp;

    }


}
