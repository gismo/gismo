/** @file gsSubdivScheme.h

    @brief Subdivision operations on mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, D. Tolis
*/


#include <gsCore/gsTemplateTools.h>

#include <gsMesh2/gsSubdivScheme.h>
#include <gsMesh2/IO.h>

#include <gsCore/gsMultiPatch.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsIO/gsXml.h>

namespace gismo {

    typedef gsSurfMesh::Point Point;
    typedef gsSurfMesh::Vertex Vertex;
    typedef gsSurfMesh::Face Face;
    typedef gsSurfMesh::Halfedge Halfedge;
    typedef gsSurfMesh::Edge Edge;

 
    void gsSubdivScheme::subdivide(index_t numSubs)
    {
        const index_t s = m_options.getInt("scheme");
        switch (s)
        {
        case 0:
            for (index_t i = 0; i != numSubs; ++i)
                cc_subdivide();
            return;
        case 1:
            for (index_t i = 0; i != numSubs; ++i)
                ds_subdivide();
            return;
        case 2:
            for (index_t i = 0; i != numSubs; ++i)
                loop_subdivide();
            return;
        }
        GISMO_ERROR("Unknown scheme.");
    }


    gsSurfMesh::Vertex_property<Point> gsSubdivScheme::vertex_limits()
    {
        const index_t s = m_options.getInt("scheme");
        switch (s)
        {
        case 0:
            return cc_vertex_limits();
        } // TODO: for the rest of subdivision schemes
        GISMO_ERROR("Unknown scheme.");
    }

    gsSurfMesh::Vertex_property<Point> gsSubdivScheme::tangent_vertex_limits()
    {
        const index_t s = m_options.getInt("scheme");
        switch (s)
        {
        case 0:
            return cc_tangent_vertex_limits();
        } // TODO: for the rest of subdivision schemes
        GISMO_ERROR("Unknown scheme.");
    }

    gsSurfMesh::Vertex_property<Point> gsSubdivScheme::normals_vertex_limits()
    {
        const index_t s = m_options.getInt("scheme");
        switch (s)
        {
        case 0:
            return cc_normals_vertex_limits();
        } // TODO: for the rest of subdivision schemes
        GISMO_ERROR("Unknown scheme.");
    }



    void gsSubdivScheme::cc_subdivide()
    {
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

    gsSurfMesh::Vertex_property<Point>
        gsSubdivScheme::cc_vertex_limits(std::string label)
    {
        auto points = m_mesh->get_vertex_property<Point>("v:point");
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


    gsSurfMesh::Vertex_property<gsSurfMesh::Point>
    gsSubdivScheme::cc_normals_vertex_limits(std::string label, bool normalize)
    {
        auto points = m_mesh->get_vertex_property<Point>("v:point");
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

    gsSurfMesh::Vertex_property<gsSurfMesh::Point>
    gsSubdivScheme::cc_tangent_vertex_limits(std::string label, bool normalize)
    {
        gsSurfMesh::Vertex v;
        gsSurfMesh::Halfedge h2;

        auto points = m_mesh->get_vertex_property<Point>("v:point");
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

    gsSurfMesh::Face_property<Point>
        gsSubdivScheme::ds_vertex_limits(std::string label)
    {
        auto limits = m_mesh->add_face_property<Point>((label), Point(0, 0, 0));
#   pragma omp parallel for default(shared)
        for (auto fit = m_mesh->faces_begin(); fit < m_mesh->faces_end(); ++fit)
        {
            if (m_mesh->is_boundary(*fit)) // TODO: Deprecate this feature by using middle edge point
            {
                gsWarn << "Boundary face is ignored.\n";
                continue;
            }
            limits[*fit] = m_mesh->face_barycenter(*fit);
        }

        return limits;
    }

    gsSurfMesh::Face_property<Point>
        gsSubdivScheme::ds_normals_vertex_limits(std::string label,
            bool normalize)
    {
        auto points = m_mesh->get_vertex_property<Point>("v:point");
        auto limits = m_mesh->add_face_property<Point>((label), Point(0, 0, 0));
        unsigned int n;
        Halfedge he,hh;
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
    void
    gsSubdivScheme::
    gradBoundary()
    {

        // Current implementation only for regular boundary (vertex valence = 3) and 
        // conrners (vertex valence = 2).
        // TODO: General case for EF in boundary by using Chebysev points (see A.Nashri 1987).

        std::map<Vertex, Point> bvmap; // New positions for boundary vertices
        Vertex bv;
        auto pts = m_mesh->points();
        // Compute the new positions for boundary vertices.
        for (auto hit : m_mesh->halfedges())
        {
            if (m_mesh->touches_boundary(hit))
            {
                bv = m_mesh->from_vertex(hit);

                if (m_mesh->valence(bv) == 3) // Regular boundary case
                {
                    bvmap[bv] = 2 * pts[bv] - pts[m_mesh->from_vertex(m_mesh->prev_halfedge(hit))];
                }
                else if (m_mesh->valence(bv) == 2) // Corner boundary case
                {
                    bvmap[bv] = 4 * pts[bv] - 2 * pts[m_mesh->from_vertex(m_mesh->prev_halfedge(hit))]
                        - 2 * pts[m_mesh->to_vertex(hit)] + pts[m_mesh->to_vertex(m_mesh->next_halfedge(hit))];
                }
                else // irregular case
                {
                    gsWarn << "Irregular boundary stop process\n";
                    return;
                }

            }
        }

        // Modify mesh boundary
        for (auto vit : m_mesh->vertices())
            if (m_mesh->is_boundary(vit))
                m_mesh->position(vit) = bvmap[vit];


    }
    void gsSubdivScheme::chaikin_scheme()
    {

        // Check if we have only curve (one face)
        GISMO_ASSERT(m_mesh->n_faces() == 1, "For a curve subdivision scheme we need one face");

        // Splitting phase (adding new vertices at edges midpoints)
        gsSurfMesh::Point tmp, tmps, tmpe;
        gsSurfMesh::Halfedge he;
        for (auto eit : m_mesh->edges())
        {
            he = m_mesh->halfedge(eit, 0);
            tmp = 0.5 * (m_mesh->position(m_mesh->from_vertex(he)) + m_mesh->position(m_mesh->to_vertex(he)));
            m_mesh->insert_vertex(eit, tmp);
        }

        he = m_mesh->halfedge(gsSurfMesh::Face(0)); // Face's halfedge

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



    void gsSubdivScheme::ds_subdivide()
    {
       
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
                    v = new_mesh.add_vertex(ds_image_point_calc_vanila(oldv, oldf));
                else
                    v = new_mesh.add_vertex(ds_image_point_calc_interpolation(oldv, oldf));

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
    gsSubdivScheme::ds_image_point_calc_interpolation(Vertex oldv, Face oldf)
    {
        unsigned int face_valence{ m_mesh->valence(oldf) };

        // Find the halfedge of the vertex I am looking in case the IDs of mesh
        // are not sequencial (i.e. Quad ID: 1,23,3,5)
        Halfedge hf = m_mesh->halfedge(oldf);
        for (auto hh : m_mesh->halfedges(oldf))
        {
            if (m_mesh->from_vertex(hh) == oldv)
            {
                hf = hh;
                break;
            }
        }

        real_t val{ 0 };

        gsSurfMesh::Point coords;
        coords.setZero();


        if (m_mesh->is_boundary(oldv)) {

            if (m_mesh->valence(oldv) == 2) { // Corner's case (same position because it will go to limit)
                coords = m_mesh->position(oldv);
            }
            else { // Chaikin method for the boundary
                if (m_mesh->is_boundary(m_mesh->to_vertex(hf))) {
                    coords = 0.75 * m_mesh->position(oldv) + 0.25 * m_mesh->position(m_mesh->to_vertex(hf));
                }
                else {
                    coords = 0.25 * m_mesh->position(m_mesh->from_vertex(m_mesh->prev_halfedge(hf))) + 0.75 * m_mesh->position(oldv);
                }
            }


        }
        else {

            int tempj{ 0 };

            // Coefficient of the first (current) vertex (i=j case)
            val = (real_t)(face_valence + 5) / (4 * face_valence);
            coords += val * m_mesh->position(m_mesh->from_vertex(hf));
            real_t sum_val{ val };
            Halfedge next_he{ m_mesh->next_halfedge(hf) };

            // Creating the mask (image vertice coefficients) by looping to the number
            // of halfedges until I reach the initial half-edge (case i!=j).
            while (next_he != hf) {
                tempj++;
                val = (3 + 2 * cos(2.0 * EIGEN_PI * (0 - tempj) / face_valence)) / (4 * face_valence);
                coords += val * m_mesh->position(m_mesh->from_vertex(next_he));
                next_he = m_mesh->next_halfedge(next_he);
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
    gsSubdivScheme::ds_image_point_calc_vanila(Vertex oldv, Face oldf)
    {
        unsigned int face_valence{ m_mesh->valence(oldf) };

        // Find the halfedge of the vertex I am looking in case the IDs of mesh
        // are not sequencial (i.e. Quad ID: 1,23,3,5)
        Halfedge hf = m_mesh->halfedge(oldf);
        for (auto hh : m_mesh->halfedges(oldf))
        {
            if (m_mesh->from_vertex(hh) == oldv)
            {
                hf = hh;
                break;
            }
        }

        real_t val{ 0 };

        gsSurfMesh::Point coords;
        coords.setZero();

        int tempj{ 0 };

        // Coefficient of the first (current) vertex (i=j case)
        val = (real_t)(face_valence + 5) / (4 * face_valence);
        coords += val * m_mesh->position(m_mesh->from_vertex(hf));
        real_t sum_val{ val };
        Halfedge next_he{ m_mesh->next_halfedge(hf) };

        // Creating the mask (image vertice coefficients) by looping to the number
        // of halfedges until I reach the initial half-edge (case i!=j).
        while (next_he != hf) {
            tempj++;
            val = (3 + 2 * cos(2.0 * EIGEN_PI * (0 - tempj) / face_valence)) / (4 * face_valence);
            coords += val * m_mesh->position(m_mesh->from_vertex(next_he));
            next_he = m_mesh->next_halfedge(next_he);
            sum_val += val;

        }

        Point temp;
        temp[0] = coords(0);
        temp[1] = coords(1);
        temp[2] = coords(2);

        return temp;

    }

    void gsSubdivScheme::loop_subdivide()
    {

        index_t mask_option = m_options.askInt("loop.maskType");

        gsSurfMesh nm;
        Vertex v;

        // Tessalate the current mesh, in order to have pure triangle mesh
        for (auto fit : m_mesh->faces())
            if (m_mesh->valence(fit) > 3)
                m_mesh->triangulate(fit);


        // Export the triangle mesh
        m_mesh->write("mesh_in_triang.off");
        // reserve vertices, edges, faces
        m_mesh->reserve(m_mesh->n_vertices() + m_mesh->n_edges() + m_mesh->n_faces(),
            2 * m_mesh->n_edges(), 4 * m_mesh->n_faces());


        index_t fnv = m_mesh->n_vertices(); // Number of vertices in the current mesh
        index_t num_f = m_mesh->n_faces(); // Number of faces in the current mesh

        gsSurfMesh::Halfedge he,nh, ph, hh;
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

        for (auto he : m_mesh->halfedges()) 
        {

            if (m_mesh->is_boundary(he)) 
            {
                continue;
            }

            v = m_mesh->from_vertex(he);
            ph = m_mesh->prev_halfedge(he);
            nm.add_triangle(innerVerts[v], edgeVerts[m_mesh->edge(he)], edgeVerts[m_mesh->edge(ph)]);
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
