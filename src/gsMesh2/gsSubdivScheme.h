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

    gsOptionList m_options;


public:

    /// Constructor accepting a reference to a mesh
    explicit gsSubdivScheme(gsSurfMesh& mesh) : m_mesh(&mesh)
    {
        m_options.addInt("scheme", "0: CC, 1: DS, 2:  Loop", 0);
        m_options.addInt("ds_opt", "Option for masks in Doo-Sabin subdivision scheme", 0);
        m_options.addInt("loop_opt", "Option for masks in Loop subdvision scheme", 0);
    }

    /// Returns the options
    gsOptionList& options() { return m_options; }

    /// Returns the mesh
    gsSurfMesh & mesh() { return *m_mesh; }

    /// Returns the mesh
    const gsSurfMesh & mesh() const { return *m_mesh; }

    void subdivide(index_t numSubs = 0)
    {
        const index_t s = m_options.getInt("scheme");
        switch(s)
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

    //todo: void vertex_limits()
    //todo: void vertex_tangents()

public: // Catmull-Clark functions
    
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

public: // Doo-Sabin functions

    /** Doo-Sabin subdvision
     * Options:
     *
     * ds_opt:
     *   0 - Interpolation of the boundary using Chaikin's scheme.
     *   1 - Vanila version that leads to trimmed boundaries.
     *
    */
    void ds_subdivide();

private:
    
    /// Doo-Sabin Image point calculation per vertex in a face (boundary interpolation)
    Point ds_image_point_calc_interpolation(Vertex oldv, Face oldf);

    /// Doo-Sabin Image point calculation per vertex in a face (trimmed)
    Point ds_image_point_calc_vanila(Vertex oldv, Face oldf);

public: // Loop subdivision

    /** Loop subdvision
     * Options:
     *  +++ 
     *
    */
    void loop_subdivide();
    


};//namespace internal

} // namespace gismo
