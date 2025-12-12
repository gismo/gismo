/** @file gsSubdivScheme.h

    @brief Subdivision operations on mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, D. Tolis
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsIO/gsXml.h>

//#define Eigen gsEigen
//EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(gismo::Point)
//#undef Eigen

#include <gsMesh2/gsProperty.h>
#include <gsMesh2/gsSurfMesh.h>
#include <gsIO/gsOptionList.h>

namespace gismo
{

/// A class for subdivision schemes in polygonal meshes.
class GISMO_EXPORT gsSubdivScheme
{

      // type definitions
    typedef gsSurfMesh::Point Point;
    typedef gsSurfMesh::Vertex Vertex;
    typedef gsSurfMesh::Face Face;

    gsSurfMesh* m_mesh;///<pointer to the input mesh

    gsOptionList s_options;


public: // Catmull-Clark functions

    explicit gsSubdivScheme(gsSurfMesh& mesh) : m_mesh(&mesh)
    { 
    
       s_options.addInt("ds_opt", "Option for masks in Doo-Sabin subdivision scheme", 0);
       s_options.addInt("loop_opt", "Option for masks in Loop subdvision scheme", 0);
    }

    gsOptionList& options() { return s_options; }
    
    gsSurfMesh mesh() { return *m_mesh; }


    /// Catmull-Clark subdivision
    void cc_subdivide();

    /// Compute CC vertex limit positions
    gsSurfMesh::Vertex_property<Point> cc_limit_points(std::string label = "v:limit");

    /// Compute CC vertex limit normals
    gsSurfMesh::Vertex_property<Point> cc_limit_normals(std::string label = "v:normal",
                                            bool normalize = true);

    /// Compute CC vertex limit tangent
    gsSurfMesh::Vertex_property<Point> cc_limit_tangent_vec(std::string label = "v:tanvec",
                                                bool normalize = true);

    // Returns true if there is a halfedge with hflag set to true emenating from vertex \a v
    inline bool has_flag(Vertex v, const gsSurfMesh::Halfedge_property<bool> & hflag);



public: // Doo-Sabin functions

    /// Doo-Sabin subdivision
    /* Doo-Sabin subdvision
     * Options:
     *
     * \t 0 - interpolation in boundary using Chaikin's scheme.
     * \t 1 - vanila version that leads to trimmed boundaries.
     *
    */
    void ds_subdivide();

    /// Doo-Sabin Image point caluculation per vertex in a face (boundary interpolation)
    Point ds_image_point_calc_interpolation(Vertex oldv, Face oldf);

    /// Doo-Sabin Image point caluculation per vertex in a face (trimmed)
    Point ds_image_point_calc_vanila(Vertex oldv, Face oldf);

public: // Loop subdivision

    ///  Loop subdivision
    /** Loop subdvision
     * Options:
     *
     * \t 0 - Simplified Loop's scheme.
     * \t 1 - Original Loop's scheme.
     *
    */
    void loop_subdivide();
    


};//namespace internal

} // namespace gismo
