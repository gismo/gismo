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

    template<size_t N>
    gsFreeformFaceData<N>::gsFreeformFaceData(const gsSurfMesh& mesh, gsSurfMesh::Face face)
    : control_points(),
      face(face)
    {
        // Create a vector of the 4 corners.
        std::vector<gismo::gsVector3d<>> corners;
        corners.reserve(4);
        for (auto const & v : mesh.vertices(face)){
            corners.emplace_back(mesh.position(v));
        }
        assert(points.size() == 4);

        // Choose the control points (N*N total) as appropriate linear combinations of the corners.
        for(size_t i = 0; i<N; i++){
          for(size_t j = 0; j<N; j++){
            real_t denom = real_t((N-1)*(N-1));
            real_t n_r = real_t(N);
            real_t i_r = real_t(i);
            real_t j_r = real_t(j);
            this->control_points(i,j) = 
                corners[0] * ((n_r-1)-j_r) * ((n_r-1)-i_r)/denom +
                corners[1] * j_r * ((n_r-1)-i_r)/denom +
                corners[3] * ((n_r-1)-j_r) * i_r/denom +
                corners[2] * j_r * i_r/denom;
          }
        }
      
    }

    template<size_t N>
    gsVector3d<>* gsFreeformFaceData<N>::vertex_control_point(
      gsSurfMesh& mesh,
      Vertex v,
      size_t inset
    ) 
    {
      // find a halfedge on the face starting on this vertex
      Halfedge hedge;
      for(auto const & he : mesh.halfedges(face)){
          if(mesh.from_vertex(he) == v) hedge = he;
      }
      auto ecp = edge_control_points(mesh, hedge, inset);
      if(ecp.size() > 0){
        return ecp[0];
      } else {
        return nullptr;
      }
    }

    template<size_t N>
    std::vector<gsVector3d<>*> gsFreeformFaceData<N>::edge_control_points(
      gsSurfMesh& mesh,
      Halfedge hedge,
      size_t inset
    )
    {
      // find the edge on the face
      size_t hedge_index(0);
      for(auto const & he : mesh.halfedges(face)){
          if(he == hedge) break;
          ++hedge_index;
      }
      // make sure it was found (does nothing in release mode)
      assert(hedge_index<4);

      std::vector<gsVector3d<>*> result;

      // We need to collect a total of `n - 2 inset` points, skipping the first and last `inset` points.
      for(size_t i = inset; i < N - inset; i++){
        switch (hedge_index) {
          case 1:
            // Edge 3: A row on the top, left to right
            result.emplace_back(&control_points(inset,i));
            break;
          case 2:
            // Edge 2: A column on the right, top to bottom
            result.emplace_back(&control_points(i,(N-1)-inset));
            break;
          case 3:
            // Edge 3: A row on the bottom, right to left
            result.emplace_back(&control_points((N-1)-inset,(N-1)-i));
            break;
          case 0:
            // Edge 0: A column on the left, bottom to top
            result.emplace_back(&control_points((N-1)-i,inset));
            break;          
        }
      }

      return result;    
    }
 
    template<size_t N>
    const gismo::gsTensorBSpline<2, real_t> gsFreeformFaceData<N>::patch() const
    {
      // Create a spline basis for a normal bezier patch.
      gsKnotVector<> kv1(0, 1, 0, N);
      gsKnotVector<> kv2(0, 1, 0, N);
      gsTensorBSplineBasis<2, real_t> basis(kv1, kv2);
      // Create a coefficient matrix out of the control points.
      //Technically, you could use just [i] and one loop here since the elements of a matrix are layed out row-wise, but this might be clearer to read.
      gsMatrix<> coeffs(N * N,3);
      for(size_t i = 0; i < N; ++i){
        for(size_t j = 0; j < N; ++j){
          int total_index = i * N + j;
          coeffs(total_index, 0) = control_points(i, j).x();
          coeffs(total_index, 1) = control_points(i, j).y();
          coeffs(total_index, 2) = control_points(i, j).z();
        }
      }

      return gsTensorBSpline<2>(basis, coeffs);
    }


    template<size_t N>
    std::array<gsMatrix<gsVector3d<>, Dynamic, Dynamic>, 2> gsFreeformSubdivision<N>::deCasteljau(const gsMatrix<gsVector3d<>, Dynamic, Dynamic>& control_net){

      // Create the 3d data vector and ensure it has the right size.
      std::vector<gsMatrix<gsVector3d<>, Dynamic, Dynamic>> points;
      points.resize(N);
      // The first layer is just the starting points
      points[0] = control_net;
      // each further layer is one shorter than the previous, but just as wide.
      for(size_t k = 1; k < N; ++k){
        points[k].resize(N-k, N);
      }

      // now construct each layer from the previous one by linear combination of adjacent points into the next layer.
      for(size_t k = 1; k < N; ++k){
        for(size_t i = 0; i < N-k; ++i){
          for(size_t j = 0; j < N; ++j){
            points[k](i,j) = (points[k-1](i,j) + points[k-1](i+1,j)) * 0.5;
          }
        }
      }

      // finally collect the first vertical layer and last-in-each-row diagonal layer into two result matrices.
      gsMatrix<gsVector3d<>, Dynamic, Dynamic> result1;
      gsMatrix<gsVector3d<>, Dynamic, Dynamic> result2;
      result1.resize(N,N);
      result2.resize(N,N);

      for(size_t i = 0; i<N; ++i){
        for(size_t j = 0; j<N; ++j){
          result1(i,j) = points[i](0,j);
          result2(i,j) = points[(N-1)-i](i,j);
        }
      }

      return {result1, result2};
    }
    
    template<size_t N>
    void gsFreeformSubdivision<N>::subdivide(gsSurfMesh &mesh)
    {
      // Remember the first vertex of each face (this is where the control nets of each face data are oriented on).
      std::vector<Vertex> first_vertices;
      for(Face f : mesh.faces()){
        first_vertices.emplace_back(mesh.to_vertex(mesh.halfedge(f)));
      }

      // Split each face into 4 and get info about the way they were split.
      std::map<gsSurfMesh::Face, std::vector<gsSurfMesh::Face>> face_map = mesh.quad_split();

      // Get face data
      gsProperty<gsFreeformFaceData<N>> face_data_vec(mesh.get_face_property<gsFreeformFaceData<N>>("bezier_points"));

      // Now fix the data on each face.
      for(auto const & parent_to_children_faces : face_map){
        // Get the face data and store it in a temporary dynamic 2d array.
        gsMatrix<gsVector3d<>, Dynamic, Dynamic> control_net(face_data_vec.vector()[parent_to_children_faces.first.idx()].control_points);

        // now control_net is a (n+1)*(n+1) matrix of control points (degree n)
        // Perform deCasteljau once to divide into two (n+1)*(n+1) matrices of control points.
        auto const first_split = this->deCasteljau(control_net);

        // Perform deCasteljau again on both of them, to get 4 (n+1)*(n+1) matrices of control points
        // In between, we need to transpose so we now divide in the other direction.
        auto top_split = this->deCasteljau(first_split[0].transpose());
        auto bot_split = this->deCasteljau(first_split[1].transpose());

        // re-transpose
        top_split[0] = top_split[0].transpose().eval();
        top_split[1] = top_split[1].transpose().eval();
        bot_split[0] = bot_split[0].transpose().eval();
        bot_split[1] = bot_split[1].transpose().eval();

        // rotate
        bot_split[1] = rotate(bot_split[1]);
        bot_split[0] = rotate(rotate(bot_split[0]));
        top_split[0] = rotate(rotate(rotate(top_split[0])));

        // Collate all these matrices in the correct order into an array.
        std::array<gsMatrix<gsVector3d<>, Dynamic, Dynamic>, 4> arr = {top_split[0], top_split[1], bot_split[1], bot_split[0]};

        // find the new face that contains the top_left vertex
        Vertex first_vertex = first_vertices[parent_to_children_faces.first.idx()];
        size_t first_face(0);
        for(size_t i = 0; i < 4; ++i){
          Face f(parent_to_children_faces.second[i]);
          for(auto const & v : mesh.vertices(f)){
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

        // Correct back references of face data and give them the correct control points.
        for(size_t f = 0; f < 4; ++f){
          auto data = &face_data_vec.vector()[children_faces_ordered[f].idx()];
          data->face = children_faces_ordered[f];
          data->control_points = arr[f];
        }

      }

    };

    template<size_t N>
    gsMatrix<gsVector3d<>, Dynamic, Dynamic> gsFreeformSubdivision<N>::rotate(const gsMatrix<gsVector3d<>, Dynamic, Dynamic>& mat){
      gsMatrix<gsVector3d<>, Dynamic, Dynamic> res;
      res.resize(mat.cols(), mat.rows());
      for(int i = 0; i<res.rows(); ++i){
        for(int j = 0; j<res.cols(); ++j){
          res(i,j) = mat(j, mat.cols()-1-i);
        }
      }
      return res;
    }

    template<size_t N>
    void gsFreeformSubdivision<N>::make_c1(gsSurfMesh &mesh)
    {
      //Get face data
      gsProperty<gsFreeformFaceData<N>> face_data_vec(mesh.get_face_property<gsFreeformFaceData<N>>("bezier_points"));
      // Now correct each face.
      for(Face f : mesh.faces()){
        auto * const face_data = &face_data_vec.vector()[f.idx()];

        // === INNER ===
        // stay the same

        // === EDGES ===
        // Iterator over all half edges adjacent to this one.
        for(auto const & hedge : mesh.halfedges(f)){
            // Get the two half-edges and then faces for this half-edge.
            auto hedge1 (hedge);
            auto hedge2(mesh.opposite_halfedge(hedge));
            auto face1 (mesh.face(hedge1)); // this should be the current face.
            auto face2 (mesh.face(hedge2));

            // Get the two N-2 sets of corresponding points.
            auto points1(face_data_vec.vector()[face1.idx()].edge_control_points(mesh, hedge1, 1));
            auto points2(face_data_vec.vector()[face2.idx()].edge_control_points(mesh, hedge2, 1));

            auto points_to_be_set = face_data->edge_control_points(mesh, hedge, 0);
      
            // Calculate the average of each control point with its partner on the other side and store it in the appropriate control point in the result.
            for(size_t i = 0; i<N-2; ++i){
                *points_to_be_set[i+1] = 
                *points1[i] * 0.5 + *points2[N-3-i] * 0.5;
            }
        }

        // === CORNERS ===
        // Iterate over all halfedges on this face.
        for(auto const & v_hedge : mesh.halfedges(f)){
            // Get all the vertices of this face as the starting points of such halfedges.
            auto v = mesh.from_vertex(v_hedge);

            // Prepare sum
            gsVector3d<real_t> sum(0.,0.,0.);
            real_t count(0.0);

            // Iterate over all halfedges leaving this vertex, one per face.
            for(auto const & out_hedge : mesh.halfedges(v)) {
                auto out_face = mesh.face(out_hedge);
                // For each of these other faces (represented by halfedges moving into v), sum over the closest inner control point. To find this one, we take all control points along the halfedge (with inset 1) and take the first, since this halfedge is outgoing with respect to v.
                sum += *face_data_vec.vector()[out_face.idx()].edge_control_points(mesh, out_hedge, 1)[0];
                count += 1.0;
            }
            // Store this as the corner control point of the result.
            // To find the correct control point, we take all control points along the halfedge and take the first, since the vertex in question is the from_vertex of the current halfedge.
            *face_data->edge_control_points(mesh, v_hedge, 0)[0] = sum/count;
        }
        
      } // end for over faces
    
    }

    template<size_t N>
    void gsFreeformSubdivision<N>::initialize_data(gsSurfMesh &mesh)
    {
      mesh.add_face_property(std::string("bezier_points"), gsFreeformFaceData<N>());
      gsProperty<gsFreeformFaceData<N>> patch_data = mesh.get_face_property<gsFreeformFaceData<N>>("bezier_points");
      for (auto f : mesh.faces()){
          patch_data.vector()[f.idx()] = gsFreeformFaceData<N>(mesh, f);
      }
    }

  template class gsFreeformSubdivision<5>;
  template class gsFreeformFaceData<5>;
  template class gsFreeformSubdivision<6>;
  template class gsFreeformFaceData<6>;
  template class gsFreeformSubdivision<9>;
  template class gsFreeformFaceData<9>;
    
}//namespace gismo
