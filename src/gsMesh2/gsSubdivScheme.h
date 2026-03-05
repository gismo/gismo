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

/// class for subdivision schemes in polygonal meshes.
template <class Scalar = real_t>
class GISMO_EXPORT gsSubdivScheme
{
protected:
    gsSurfMesh<Scalar>* m_mesh;///<pointer to the input mesh

    gsOptionList m_options;

public:
    // type definitions
    typedef typename gsSurfMesh<Scalar>::Point Point;
    typedef typename gsSurfMesh<Scalar>::Vertex Vertex;
    typedef typename gsSurfMesh<Scalar>::Face Face;

public:

    /// Constructor accepting a reference to a mesh
    explicit gsSubdivScheme(gsSurfMesh<Scalar>& mesh) : m_mesh(&mesh)
    {
        m_options.addInt("scheme", "0: CC, 1: DS, 2:  Loop", 0);
        m_options.addInt("ds.boundaryMask", "Option for masks in Doo-Sabin subdivision scheme",0);
        m_options.addInt("loop.maskType",   "Option for masks in Loop subdvision scheme", 0);
    }

    /// Returns the options
    gsOptionList& options() { return m_options; }

    /// Returns the mesh
    gsSurfMesh<Scalar> & mesh() { return *m_mesh; }

    /// Returns the mesh
    const gsSurfMesh<Scalar> & mesh() const { return *m_mesh; }

    /// Generic method for subdivide
    void subdivide(index_t numSubs = 0);

    /// Method for calculating vertex in limit of subdivision scheme
    typename gsSurfMesh<Scalar>::template Vertex_property<Point> vertex_limits();
    
    /// Method for calculating tangents on vertices in limit of subdivision scheme
    typename gsSurfMesh<Scalar>::template Vertex_property<Point> tangent_vertex_limits();
    
    /// Method for calculating normals on vertices in limit of subdivision scheme
    typename gsSurfMesh<Scalar>::template Vertex_property<Point> normals_vertex_limits();

public: // Catmull-Clark functions
    
    /// Catmull-Clark subdivision
    void cc_subdivide();

    /// Compute CC vertex limit positions
    typename gsSurfMesh<Scalar>::template Vertex_property<Point> cc_vertex_limits(std::string label = "v:limit");

    /// Compute CC vertex limit normals
    typename gsSurfMesh<Scalar>::template Vertex_property<Point> cc_normals_vertex_limits(std::string label = "v:normal",
                                            bool normalize = true);

    /// Compute CC vertex limit tangent
    typename gsSurfMesh<Scalar>::template Vertex_property<Point> cc_tangent_vertex_limits(std::string label = "v:tanvec",
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
     *
     * loop_opt:
     *   0 - Simplified Loop's scheme. (cf. book Warren, Weimer 2002)
     *   1 - Original Loop's scheme.  (cf. book Loop 1987)
     *
    */
    void loop_subdivide();
    


};//namespace internal

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsSubDivScheme.hpp)
#endif
