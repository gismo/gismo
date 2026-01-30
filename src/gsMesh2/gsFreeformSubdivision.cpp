/** @file gsFreeformSubdivision.cpp

    @brief Classes for Freeform Subdivision on a quadrangular mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Mussmaecher
*/

#include "gsCore/gsDebug.h"
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

    gsVector3d<>* gsFreeformFaceData::vertex_control_point(
      const gsSurfMesh& mesh,
      Vertex v,
      size_t inset
    ){
      // find a halfedge on the face starting on this vertex
      Halfedge hedge;
      for(auto he : mesh.halfedges(face)){
          if(mesh.from_vertex(he) == v) hedge = he;
      }
      return gsFreeformFaceData::edge_control_points(mesh, hedge, inset)[0];
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


    std::array<gsMatrix<gsVector3d<>, Dynamic, Dynamic>, 2> gsFreeformSubdivision::deCasteljau(const gsMatrix<gsVector3d<>, Dynamic, Dynamic>& patch_vec){

      int n = patch_vec.cols();

      // Create the 3d data vector and ensure it has the right size.
      std::vector<gsMatrix<gsVector3d<>, Dynamic, Dynamic>> points;
      points.resize(n);
      // The first layer is just the starting points
      points[0] = patch_vec;
      // each further layer is one shorter than the previous, but just as wide.
      for(int k = 1; k < n; ++k){
        points[k].resize(n-k, n);
      }

      // now construct each layer from the previous one by linear combination of adjacent points into the next layer.
      for(int k = 1; k <= n; ++k){
        for(int i = 0; i < n-k; ++i){
          for(int j = 0; j < n; ++j){
            points[k](i,j) = (points[k-1](i,j) + points[k-1](i+1,j)) * 0.5;
          }
        }
      }

      // finally collect the first vertical layer and last-in-each-row diagonal layer into two result matrices.
      gsMatrix<gsVector3d<>, Dynamic, Dynamic> result1;
      gsMatrix<gsVector3d<>, Dynamic, Dynamic> result2;
      result1.resize(n,n);
      result2.resize(n,n);

      for(int i = 0; i<n; ++i){
        for(int j = 0; j<n; ++j){
          result1(i,j) = points[i](0,j);
          result2(i,j) = points[(n-1)-i](i,j);
        }
      }

      return {result1, result2};
    }
    
    void gsFreeformSubdivision::subdivide(gsSurfMesh &mesh)
    {
      // Remember the first vertex of each face (this is where the patches are oriented on).
      std::vector<Vertex> first_vertices;
      for(Face f : mesh.faces()){
        first_vertices.emplace_back(mesh.to_vertex(mesh.halfedge(f)));
      }

      // Split each face into 4 and get info about the way they were split.
      std::map<gsSurfMesh::Face, std::vector<gsSurfMesh::Face>> face_map = mesh.quad_split();

      // Get patch data
      gsProperty<gsFreeformFaceData> patch_data(mesh.get_face_property<gsFreeformFaceData>("bezier_points"));

      // Now fix the data on each face.
      for(auto parent_to_children_faces : face_map){
        // Get the patch data and store it in a temporary dynamic 2d array.
        gsMatrix<gsVector3d<>, Dynamic, Dynamic> patch(patch_data.vector()[parent_to_children_faces.first.idx()].control_points);

        // now patch is a (n+1)*(n+1) matrix of control points (degree n)
        // Perform deCasteljau once to divide into two (n+1)*(n+1) matrices of control points.
        auto cjtemp = this->deCasteljau(patch);

        // Perform deCasteljau again on both of them, to get 4 (n+1)*(n+1) matrices of control points
        // In between, we need to transpose so we now divide in the other direction.
        auto cj12 = this->deCasteljau(cjtemp[0].transpose());
        auto cj34 = this->deCasteljau(cjtemp[1].transpose());

        // re-transpose
        cj12[0] = cj12[0].transpose().eval();
        cj12[1] = cj12[1].transpose().eval();
        cj34[0] = cj34[0].transpose().eval();
        cj34[1] = cj34[1].transpose().eval();

        // rotate
        cj34[1] = rotate(cj34[1]);
        cj34[0] = rotate(rotate(cj34[0]));
        cj12[0] = rotate(rotate(rotate(cj12[0])));

        // Collate all these matrices in the correct order into an array.
        std::array<gsMatrix<gsVector3d<>, Dynamic, Dynamic>, 4> arr = {cj12[0], cj12[1], cj34[1], cj34[0]};

        // find the new face that contains the top_left vertex
        Vertex first_vertex = first_vertices[parent_to_children_faces.first.idx()];
        size_t first_face(0);
        for(size_t i = 0; i < 4; ++i){
          Face f(parent_to_children_faces.second[i]);
          for(Vertex v : mesh.vertices(f)){
            if(v == first_vertex){
              first_face = i;
            }
          }
        }

        // Collate the faces into a correctly ordered array as well.
        std::array<Face, 4> children_faces_ordered;
        for(int i = 0; i < 4; ++i){
          children_faces_ordered[i] = parent_to_children_faces.second[(i + first_face)%4];
        }

        // Correct back references of patch data and give them the correct control points.
        for(size_t f = 0; f < 4; ++f){
          auto data = &patch_data.vector()[children_faces_ordered[f].idx()];
          data->face = children_faces_ordered[f];
          data->control_points = arr[f];
        }

      }

    };

    gsMatrix<gsVector3d<>, Dynamic, Dynamic> rotate(const gsMatrix<gsVector3d<>, Dynamic, Dynamic>& mat){
      gsMatrix<gsVector3d<>, Dynamic, Dynamic> res;
      res.resize(mat.cols(), mat.rows());
      for(int i = 0; i<res.rows(); ++i){
        for(int j = 0; j<res.cols(); ++j){
          res(i,j) = mat(j, mat.cols()-1-i);
        }
      }
      return res;
    }

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
