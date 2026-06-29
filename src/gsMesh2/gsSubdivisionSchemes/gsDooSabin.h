/** @file gsDooSabin.h

    @brief Doo-Sabin subdivision on a  mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, D. Tolis, L. Mussmaecher
*/

#pragma once

#include <gsMesh2/gsSubdivisionScheme.h>

namespace gismo
{

/// class for subdivision schemes in polygonal meshes.
class GISMO_EXPORT gsDooSabin : public gsSubdivisionScheme
{

public: // Constructors

    /// \brief Constructor with a mesh to target.
    ///
    /// Constructor that accepts a mesh to be targeted by this constructor.
    /// Creates the 'ds.boundaryMask' option and initializes it with value `0`.
    /// 
    /// ds.boundaryMask (Boundary treatment):
    ///  * 0 - interpolatory case (cf. Chaikin 1974). 
    ///  * 1 - trimmed boundary (cf. Doo-Sabin 1978). 
    explicit gsDooSabin(gsSurfMesh* mesh = nullptr) : gsSubdivisionScheme()
    {
        m_options.addInt("ds.boundaryMask", "Option for boundary masks in Doo-Sabin subdivision scheme",0);
        this->assign(mesh);
    }

public:
    void subdivide_impl() GISMO_OVERRIDE;

    /// Compute face limit positions
    gsSurfMesh::Face_property<Point> face_limits(std::string label = "f:limit");

    /// Compute face limit normals
    gsSurfMesh::Face_property<Point> face_normal_limits(std::string label = "f:normal");

    /// Compute face limit tangent
    gsSurfMesh::Face_property<Point> face_tangent_limits(std::string label = "f:tanvec")
    {GISMO_NO_IMPLEMENTATION}



private: // Helper functions

    /// Doo-Sabin Image point calculation per halfedge in a face (vanila - trimmed)
    ///
    /// \param hf: Given halfedge. Its start will be the vertex from which we will compute
    /// the image.
    Point ds_image_point_calc(Halfedge hf);

    /// Doo-Sabin Image point calculation per halfedge in a face (boundary interpolation - Chaikin scheme)
    ///
    /// \param hf: Given halfedge. Its start will be the vertex from which we will compute
    /// the image.
    Point ds_image_point_calc_interpolation(Halfedge hf);
  
};//namespace internal

} // namespace gismo

