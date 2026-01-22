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
class GISMO_EXPORT gsCatmullClark : public gsSubdivisionScheme
{

public: // Constructors
  /// Default constructor.
  /// Catmull-Clark has no special options.
  gsCatmullClark() : gsSubdivisionScheme() {}

public:
  void subdivide(gsSurfMesh *mesh) override;
  
};//namespace internal

} // namespace gismo

