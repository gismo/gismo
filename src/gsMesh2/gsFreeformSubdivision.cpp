/** @file gsFreeformSubdivision.cpp

    @brief Classes for Freeform Subdivision on a quadrangular mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Mussmaecher
*/

#include <gsNurbs/gsTensorBSpline.h>
#include <gsMesh2/gsSubdivisionScheme.h>
#include <gsMesh2/gsFreeformSubdivision.h>
#include <iterator>

namespace gismo
{

    const std::vector<std::array<size_t, 5>> gsFreeformFaceData::gs_FREEFORMDATA_EDGES = {
      {20, 15, 10, 5, 0},
      {0, 1, 2, 3, 4},
      {4, 9, 14, 19, 24},
      {24, 23, 22, 21, 20},
    };

    const std::vector<std::array<size_t, 3>> gsFreeformFaceData::gs_FREEFORMDATA_INNEREDGES = {
      {16,11,6},
      {6,7,8},
      {8,13,18},
      {18,17,16},
    };

    gsFreeformFaceData::gsFreeformFaceData(const gsSurfMesh& mesh, gsSurfMesh::Face face)
    : control_points(),
      face(face)
    {
        // Create a vector of the 4 corners.
        std::vector<gismo::gsVector3d<>> points;
        points.reserve(4);
        for (auto v : mesh.vertices(face)){
            points.emplace_back(mesh.position(v));
        }
        assert(points.size() == 4);

        // Choose the control points (out of 25 total) as appropriate linear combinations of the corners.
        for(size_t i = 0; i<25; ++i){
            real_t ih(i % 5);                // horizontal index in the grid, varies from 0-4
            real_t iv(std::floor(i / 5)); // vertical index in the grid, varies from 0-4
          
            this->control_points[i] =
                points[0] * (4-ih) * (4-iv)/16. +
                points[1] * ih * (4-iv)/16. +
                points[3] * (4-ih) * iv/16. +
                points[2] * ih * iv/16.;
        }
       
    }

    std::array<gismo::gsVector3d<>, 5> gsFreeformFaceData::edge_control_points(
      const gsSurfMesh& mesh,
      Halfedge hedge
    )
    {
      // find the edge on the face
      size_t hedge_position(0);
      for(auto he : mesh.halfedges(face)){
          if(he == hedge) break;
          ++hedge_position;
      }
      // make sure it was found (does nothing in release mode)
      assert(hedge_position<4);

      // Take the prepared indices from above to put the correct control points into a result vector.
      std::array<gismo::gsVector3d<>, 5> result;
      for(size_t i = 0;i<5;i++){
        result[i] = this->control_points[gsFreeformFaceData::gs_FREEFORMDATA_EDGES[hedge_position][i]];
      }

      return result;
    }

    std::array<gismo::gsVector3d<>, 3> gsFreeformFaceData::edge_inner_control_points(
      const gsSurfMesh& mesh,
      Halfedge hedge
    )
    {
      // find the edge on the face
      size_t hedge_position(0);
      for(auto he : mesh.halfedges(face)){
          if(he == hedge) break;
          ++hedge_position;
      }
      // make sure it was found (does nothing in release mode)
      assert(hedge_position<4);

      // Take the prepared indices from above to put the correct control points into a result vector.
      std::array<gismo::gsVector3d<>, 3> result;
      for(size_t i = 0;i<3;i++){
        result[i] = this->control_points[gsFreeformFaceData::gs_FREEFORMDATA_INNEREDGES[hedge_position][i]];
      }

      return result;
    }

    gismo::gsVector3d<> gsFreeformFaceData::vertex_control_point(
      const gsSurfMesh& mesh,
      Vertex v
    )
    {
      // find the vertex on the face
      size_t v_position(0);
      for(auto vertex: mesh.vertices(face)){
          if(vertex == v) break;
          ++v_position;
      }
      // make sure it was found (does nothing in release mode)
      assert(v_position<4);

      // Return the correct result based on the prepared index above.

      return control_points[gsFreeformFaceData::gs_FREEFORMDATA_EDGES[v_position][4]];
    }

    gismo::gsVector3d<> gsFreeformFaceData::vertex_inner_control_point(
      const gsSurfMesh& mesh,
      Vertex v
    )
    {
      // find the vertex on the face
      size_t v_position(0);
      for(auto vertex: mesh.vertices(face)){
          if(vertex == v) break;
          ++v_position;
      }
      // make sure it was found (does nothing in release mode)
      assert(v_position<4);

      // Return the correct result based on the prepared index above.

      return control_points[gsFreeformFaceData::gs_FREEFORMDATA_INNEREDGES[v_position][2]];
    }

    gismo::gsTensorBSpline<2, real_t> gsFreeformFaceData::patch()
    {
      // Create a spline basis for a normal bezier patch.
      gsKnotVector<> kv1(0, 1, 0, 5);
      gsKnotVector<> kv2(0, 1, 0, 5);
      gsTensorBSplineBasis<2, real_t> basis(kv1, kv2);
      // Create a coefficient matrix out of the control points.
      gsMatrix<> coeffs(25,3);
      for(size_t i = 0; i < control_points.size(); ++i){
          coeffs(i, 0) = control_points[i].x();
          coeffs(i, 1) = control_points[i].y();
          coeffs(i, 2) = control_points[i].z();
      }

      return gsTensorBSpline<2>(basis, coeffs);
    }


    
    void gsFreeformSubdivision::subdivide(gsSurfMesh *mesh)
    {
      //TODO
    };

    void gsFreeformSubdivision::make_c1(gsSurfMesh &mesh)
    {
      //Get Patch data
      gsProperty<gsFreeformFaceData> patch_data = mesh.get_face_property<gsFreeformFaceData>("bezier_points");
      // Now correct each face.
      for(Face f : mesh.faces()){
        auto patch = &patch_data.vector()[f.idx()];

        // === INNER ===
        // stay the same

        // === EDGES ===
        // Iterator over all half edges adjacent to this one.
        size_t hedge_counter(0);
        for(auto hedge : mesh.halfedges(f)){
            // Get the two half-edges and then faces for this half-edge.
            auto hedge1 (hedge);
            auto hedge2(mesh.opposite_halfedge(hedge));
            auto face1 (mesh.face(hedge1)); // this should be the current face.
            auto face2 (mesh.face(hedge2));

            // Get the two 3-sets of corresponding points.
            auto points1 = patch_data.vector()[face1.idx()].edge_inner_control_points(mesh, hedge1);
            auto points2 = patch_data.vector()[face2.idx()].edge_inner_control_points(mesh, hedge2);
          
            // Calculate the average of each control point with its partner on the other side and store it in the appropriate control point in the result.
            for(int i = 0; i<3; ++i){
                patch->control_points[gsFreeformFaceData::gs_FREEFORMDATA_EDGES[hedge_counter][i+1]] = points1[i] * 0.5 + points2[2-i] * 0.5;
            }

            ++hedge_counter;
        }

        // === CORNERS ===
        // Prepare indices.
        std::vector<size_t> vertex_result = {
            0, 4, 24, 20
        };

        // Iterate over all faces adjacent to this vertex.
        size_t vertex_counter(0);
        for(auto v : mesh.vertices(f)){
            // Sum over the inner control point neares to this corner for each face.
            gsVector3d<real_t> sum(gismo::gsVector3d<real_t>(0., 0., 0.));
            real_t count(0.0);
            for(auto extra_face : mesh.faces(v)) {
                sum += patch_data.vector()[extra_face.idx()].vertex_inner_control_point(mesh, v);
                count += 1.0;
            }
            // Store this as the corner control point of the result.
            patch->control_points[gsFreeformFaceData::gs_FREEFORMDATA_EDGES[vertex_counter][4]] = sum / count;
            ++vertex_counter;
        }
        
      } // end for over faces
    
    }
    
}//namespace gismo
