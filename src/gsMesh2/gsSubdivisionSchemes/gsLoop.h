/** @file gsLoop.h

    @brief Loop subdivision on a triangular  mesh.

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
class GISMO_EXPORT gsLoop : public gsSubdivisionScheme
{

public: // Constructors
  /// \brief Default constructor.
  ///
  /// Default constructor.
  /// Creates the 'loop.maskType' optionw and initializes it with value `0`.
  gsLoop() : gsSubdivisionScheme()
  {
    m_options.addSwitch("normalize", "Normalize limit normals and tangents", true);  
    m_options.addInt("loop.maskType", "Option for mask in Loop subdivision scheme",1);
  }

  /// \brief Constructor with a mesh to target.
  ///
  /// Constructor that accepts a mesh to be targeted by this constructor.
  /// Creates the 'loop.maskType' option and initializes it with value `1`.
  /// 
  /// loop.maskType:
  ///  * 0 - Simplified Loop's scheme. (cf. book Warren, Weimer 2002) 
  ///  * 1 - Original Loop's scheme.  (cf. book Loop 1987)
  gsLoop(gsSurfMesh* mesh) : gsSubdivisionScheme()
  {
      m_options.addInt("loop.maskType", "Option for mask in Loop subdivision scheme",1);
      m_options.addSwitch("normalize", "Normalize limit normals and tangents", true);
      this->assign(mesh);
  }

public:
  void subdivide_impl() GISMO_OVERRIDE;

private: // Helper functions

    
  
};//namespace internal

} // namespace gismo

