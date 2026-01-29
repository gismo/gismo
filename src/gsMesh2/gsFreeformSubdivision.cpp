/** @file gsFreeformSubdivision.cpp

    @brief Classes for Freeform Subdivision on a quadrangular mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Mussmaecher
*/

#include "gsCore/gsDebug.h"
#include <cstddef>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsMesh2/gsSubdivisionScheme.h>
#include <gsMesh2/gsFreeformSubdivision.h>

namespace gismo
{

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

        // Choose the control points (25 total) as appropriate linear combinations of the corners.
        for(int i = 0; i<control_points.rows(); i++){
          for(int j = 0; j<control_points.cols(); j++){
            this->control_points(i,j) = 
                points[0] * (4-j) * (4-i)/16. +
                points[1] * j * (4-i)/16. +
                points[3] * (4-j) * i/16. +
                points[2] * j * i/16.;
          }
        }
      
    }

    std::vector<gsVector3d<>*> gsFreeformFaceData::edge_control_points(
      const gsSurfMesh& mesh,
      Halfedge hedge,
      size_t inset
    ){
      // find the edge on the face
      size_t hedge_position(0);
      for(auto he : mesh.halfedges(face)){
          if(he == hedge) break;
          ++hedge_position;
      }
      // make sure it was found (does nothing in release mode)
      assert(hedge_position<4);

      size_t n = control_points.cols();
      std::vector<gsVector3d<>*> result;

      // We need to collect a total of `n - 2 inset` points, skipping the first and last `inset` points.
      for(size_t i = inset; i < n - inset; i++){
        switch (hedge_position) {
          case 1:
            // Edge 3: A row on the top, left to right
            result.emplace_back(&control_points(inset,i));
            break;
          case 2:
            // Edge 2: A column on the right, top to bottom
            result.emplace_back(&control_points(i,(n-1)-inset));
            break;
          case 3:
            // Edge 3: A row on the bottom, right to left
            result.emplace_back(&control_points((n-1)-inset,(n-1)-i));
            break;
          case 0:
            // Edge 0: A column on the left, bottom to top
            result.emplace_back(&control_points((n-1)-i,inset));
            break;          
        }
      }

      return result;    
    }
 
    gismo::gsTensorBSpline<2, real_t> gsFreeformFaceData::patch()
    {
      // Create a spline basis for a normal bezier patch.
      gsKnotVector<> kv1(0, 1, 0, 5);
      gsKnotVector<> kv2(0, 1, 0, 5);
      gsTensorBSplineBasis<2, real_t> basis(kv1, kv2);
      // Create a coefficient matrix out of the control points.
      gsMatrix<> coeffs(25,3);
      for(int i = 0; i < control_points.size(); ++i){
          //Technically, you could use just [i] here since the elements of a matrix are layed out row-wise, but this might be clearer to read.
          size_t i_row = i % 5;
          size_t i_col = i / 5;
          coeffs(i, 0) = control_points(i_row, i_col).x();
          coeffs(i, 1) = control_points(i_row, i_col).y();
          coeffs(i, 2) = control_points(i_row, i_col).z();
      }

      return gsTensorBSpline<2>(basis, coeffs);
    }


    
    void gsFreeformSubdivision::subdivide(gsSurfMesh &mesh)
    {
      // Split each face into 4 and get info about the way they were split.
      std::map<gsSurfMesh::Face, std::vector<gsSurfMesh::Face>> x = mesh.quad_split();

      // Get patch data
      gsProperty<gsFreeformFaceData> patch_data(mesh.get_face_property<gsFreeformFaceData>("bezier_points"));

      // size_t degree(4);
      // size_t degree_target(degree * 2);

      for(auto e : x){
        // print some info
        gsInfo << e.first << " -> ";
        for(auto f : e.second){
          gsInfo << f << ", ";
        }
        gsInfo << "\n";

        // Correct back references of patch data
        for(auto f : e.second){
          patch_data.vector()[f.idx()].face = f;
        }


        //TODO: insert correct patch data

      }
    };

    void gsFreeformSubdivision::make_c1(gsSurfMesh &mesh)
    {
      //Get Patch data
      gsProperty<gsFreeformFaceData> patch_data(mesh.get_face_property<gsFreeformFaceData>("bezier_points"));
      // Now correct each face.
      for(Face f : mesh.faces()){
        auto patch = &patch_data.vector()[f.idx()];

        // === INNER ===
        // stay the same

        // === EDGES ===
        // Iterator over all half edges adjacent to this one.
        for(auto hedge : mesh.halfedges(f)){
            // Get the two half-edges and then faces for this half-edge.
            auto hedge1 (hedge);
            auto hedge2(mesh.opposite_halfedge(hedge));
            auto face1 (mesh.face(hedge1)); // this should be the current face.
            auto face2 (mesh.face(hedge2));

            // Get the two 3-sets of corresponding points.
            auto points1(patch_data.vector()[face1.idx()].edge_control_points(mesh, hedge1, 1));
            auto points2(patch_data.vector()[face2.idx()].edge_control_points(mesh, hedge2, 1));
          
            // Calculate the average of each control point with its partner on the other side and store it in the appropriate control point in the result.
            for(int i = 0; i<3; ++i){
                *patch->edge_control_points(mesh, hedge, 0)[i+1] = 
                *points1[i] * 0.5 + *points2[2-i] * 0.5;
            }
        }

        // === CORNERS ===
        // Iterate over all halfedges on this face.
        for(auto v_hedge : mesh.halfedges(f)){
            // Get all the vertices of this face as the starting points of such halfedges.
            auto v = mesh.from_vertex(v_hedge);

            // Prepare sum
            gsVector3d<real_t> sum(0.,0.,0.);
            real_t count(0.0);

            // Iterate over all halfedges leaving this vertex, one per face.
            for(auto out_hedge : mesh.halfedges(v)) {
                auto out_face = mesh.face(out_hedge);
                // For each of these other faces (represented by halfedges moving into v), sum over the closest inner control point. To find this one, we take all control points along the halfedge (with inset 1) and take the first, since this halfedge is outgoing with respect to v.
                sum += *patch_data.vector()[out_face.idx()].edge_control_points(mesh, out_hedge, 1)[0];
                count += 1.0;
            }
            // Store this as the corner control point of the result.
            // To find the correct control point, we take all control points along the halfedge and take the first, since the vertex in question is the from_vertex of the current halfedge.
            *patch->edge_control_points(mesh, v_hedge, 0)[0] = sum/count;
                gsInfo<<"all down here\n";
        }
        
      } // end for over faces
    
    }
    
}//namespace gismo
