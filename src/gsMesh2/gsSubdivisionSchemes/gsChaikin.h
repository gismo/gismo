/** @file gsChaikin.h

    @brief Chaikin subdivision on a polyline.

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
class GISMO_EXPORT gsChaikin : public gsSubdivisionScheme<Scalar>
{

public: // Constructors
  /// \brief Default constructor.
  ///
  /// Default constructor.
  gsChaikin() : gsSubdivisionScheme<Scalar>(){}

  /// \brief Constructor with a mesh to target.
  ///
  /// Constructor that accepts a mesh to be targeted by this constructor.
  explicit gsChaikin(gsSurfMesh<Scalar>* mesh) : gsSubdivisionScheme<Scalar>()
  {
      this->assign(mesh);
      GISMO_ASSERT(this->m_mesh->n_faces() <= 1, "For a curve subdivision scheme we need at most one face");
  }

public:
  void subdivide_impl() GISMO_OVERRIDE;

};//namespace internal

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsChaikin.hpp)
#endif
