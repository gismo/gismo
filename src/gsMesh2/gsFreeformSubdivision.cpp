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

        size_t n = control_points.cols();

        // Choose the control points (n*n=25 total) as appropriate linear combinations of the corners.
        for(int i = 0; i<control_points.rows(); i++){
          for(int j = 0; j<control_points.cols(); j++){
            this->control_points(i,j) = 
                points[0] * ((n-1)-j) * ((n-1)-i)/(n-1)/(n-1) +
                points[1] * j * ((n-1)-i)/(n-1)/(n-1) +
                points[3] * ((n-1)-j) * i/(n-1)/(n-1) +
                points[2] * j * i/(n-1)/(n-1);
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
      size_t n = control_points.cols();
      // Create a spline basis for a normal bezier patch.
      gsKnotVector<> kv1(0, 1, 0, n);
      gsKnotVector<> kv2(0, 1, 0, n);
      gsTensorBSplineBasis<2, real_t> basis(kv1, kv2);
      // Create a coefficient matrix out of the control points.
      gsMatrix<> coeffs(n * n,3);
      for(int i = 0; i < control_points.size(); ++i){
          //Technically, you could use just [i] here since the elements of a matrix are layed out row-wise, but this might be clearer to read.
          size_t i_row = i % n;
          size_t i_col = i / n;
          coeffs(i, 0) = control_points(i_row, i_col).x();
          coeffs(i, 1) = control_points(i_row, i_col).y();
          coeffs(i, 2) = control_points(i_row, i_col).z();
      }

      return gsTensorBSpline<2>(basis, coeffs);
    }


    
    void gsFreeformSubdivision::subdivide(gsSurfMesh &mesh)
    {
      // Split each face into 4 and get info about the way they were split.
      std::map<gsSurfMesh::Face, std::vector<gsSurfMesh::Face>> face_map = mesh.quad_split();

      // Get patch data
      gsProperty<gsFreeformFaceData> patch_data(mesh.get_face_property<gsFreeformFaceData>("bezier_points"));

      // degree
      size_t n = patch_data[0].control_points.cols() - 1;

      for(auto e : face_map){

        // DeCasteljau
        auto patch = patch_data.vector()[e.first.idx()];

        // Create two temporary holders for deCasteljau
        std::vector<std::vector<gsVector3d<>>> patch_vec;
        std::vector<std::vector<gsVector3d<>>> patch_vec_2;

        // copy patch data
        for(int i = 0; i < patch.control_points.rows(); ++i){
          patch_vec.emplace_back();
          for(int j = 0; j < patch.control_points.cols(); ++j){
            patch_vec[i].emplace_back(patch.control_points(i,j));
          }          
        }

        // de Casteljau (vertical, n times)
        for(size_t k = 0; k<n; ++k){
          // 0th row
          patch_vec_2.push_back(patch_vec[0]);

          // all the other rows
          for(size_t i = 0; i+1 < patch_vec.size(); ++i){
            patch_vec_2.emplace_back();
            // for each column
            for(size_t j = 0; j < patch_vec[i].size(); ++j){
              patch_vec_2[i+1].emplace_back(patch_vec[i][j] * 0.5 + patch_vec[i+1][j] * 0.5);
            }          
          }

          // last row
          patch_vec_2.push_back(patch_vec[patch_vec.size()-1]);

          // put the new vector into the old one, then delete the new one
          patch_vec = patch_vec_2;
          patch_vec_2 = std::vector<std::vector<gsVector3d<>>>();
        }


        // de Casteljau (horizontal, n times)
        for(size_t k = 0; k<n; ++k){

          // for each row
          for(size_t i = 0; i < patch_vec.size(); ++i){
            patch_vec_2.emplace_back();
            // first column
            patch_vec_2[i].emplace_back(patch_vec[i][0]);
            // middle columns
            for(size_t j = 0; j+1 < patch_vec[i].size(); ++j){
              patch_vec_2[i].emplace_back(patch_vec[i][j] * 0.5 + patch_vec[i][j+1] * 0.5);
            }          
            // last column
            patch_vec_2[i].emplace_back(patch_vec[i][patch_vec[i].size() - 1]);
          }

          // put the new vector into the old one, then delete the new one
          patch_vec = patch_vec_2;
          patch_vec_2 = std::vector<std::vector<gsVector3d<>>>();
        }

        // now patch_vec is a (2n+1)*(2n+1) matrix of control points

        // Correct back references of patch data
        size_t face_counter(0);
        for(auto f : e.second){
          auto data = &patch_data.vector()[f.idx()];
          data->face = f;
          for(int i =0; i<data->control_points.rows(); ++i){
            for(int j =0; j<data->control_points.cols(); ++j){
              data->control_points(i,j) = patch_vec[
                i + n * (face_counter % 2)
              ][
                j + n * (face_counter / 2) 
              ];
            }
          }

          ++face_counter;
        }

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
        }
        
      } // end for over faces
    
    }
    
}//namespace gismo
