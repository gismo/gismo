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
class GISMO_EXPORT gsDooSabin : public gsSubdivisionScheme
{

public: // Constructors
  /// Default constructor.
  /// Catmull-Clark has no special options.
  gsDooSabin() : gsSubdivisionScheme() {
    m_options.addInt("ds.boundaryMask", "Option for masks in Doo-Sabin subdivision scheme",0);
  }

public:
  void subdivide(gsSurfMesh *mesh) override;

private: // Helper functions
    
    /// Doo-Sabin Image point calculation per vertex in a face (boundary interpolation)
    Point ds_image_point_calc_interpolation(Vertex oldv, Face oldf, const gsSurfMesh& mesh);

    /// Doo-Sabin Image point calculation per vertex in a face (trimmed)
    Point ds_image_point_calc_vanila(Vertex oldv, Face oldf, const gsSurfMesh& mesh);
  
};//namespace internal

} // namespace gismo

