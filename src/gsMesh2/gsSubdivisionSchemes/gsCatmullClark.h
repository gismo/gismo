/** @file gsCatmullClark.h

    @brief Catmull-Clark subdivision on a  mesh.

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
template <class Scalar=real_t>
class GISMO_EXPORT gsCatmullClark : public gsSubdivisionScheme<Scalar>
{
public:
    typedef gsSubdivisionScheme<Scalar> Base;
    typedef typename gsSurfMesh<Scalar>::Point Point;
    typedef typename gsSurfMesh<Scalar>::Halfedge Halfedge;

public: // Constructors
    /// \brief Constructor with a mesh to target.
    ///
    /// Constructor that accepts a mesh to be targeted by this constructor.
    /// Catmull-Clark has no special options.
    explicit gsCatmullClark(gsSurfMesh<Scalar>* mesh = nullptr): gsSubdivisionScheme<Scalar>()
    {
        this->m_options.addSwitch("normalize", "Normalize limit normals and tangents", true);
        this->assign(mesh);
    }

    static void apply(gsSurfMesh<Scalar>& mesh);

    /// Compute vertex limit positions
    typename gsSurfMesh<Scalar>::template Vertex_property<Point> vertex_limits(std::string label = "v:limit");

    /// Compute vertex limit normals for Catmull-Clark subdivision scheme
    typename gsSurfMesh<Scalar>::template Vertex_property<Point> vertex_normal_limits(std::string label = "v:normal");

    /// Compute vertex limit tangent for Catmull-Clark subdivision scheme
    typename gsSurfMesh<Scalar>::template Vertex_property<Point> vertex_tangent_limits(std::string label = "v:tanvec");

protected:

    void subdivide_impl() GISMO_OVERRIDE;

};//namespace internal

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsCatmullClark.hpp)
#endif
