/** @file gsReadOpenNurbs.h

    @brief Declaration of function for data input from the Rhinoceros 3DM file format.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsDebug.h>

#include <gsIO/gsXml.h>

class ON_Surface;
class ON_Curve;
class ON_NurbsSurface;
class ON_NurbsCurve;
class ON_NurbsCage;
class ON_Brep;
class ON_MorphControl;

namespace gismo {

namespace extensions {

  /// Reads a 3dm file using OpenNurbs
  bool gsReadOpenNurbs( const char * arg, internal::gsXmlTree & data  );
  
  /// Extracts a surface from OpenNurbs
  bool readON_Surface     ( const ON_Surface * pusrface, internal::gsXmlTree & data  );

  /// Extracts a curve from OpenNurbs
  bool readON_Curve       ( const ON_Curve   * pcurve  , internal::gsXmlTree & data  );

  /// Extracts a NURBS surface from OpenNurbs
  bool readON_NurbsSurface( const ON_NurbsSurface * psurface, internal::gsXmlTree & data  );

  /// Extracts a NURBS curve from OpenNurbs
  bool readON_NurbsCurve  ( const ON_NurbsCurve   * pcurve  , internal::gsXmlTree & data  );

  /// Extracts a NURBS cage from OpenNurbs
  bool readON_NurbsCage   ( const ON_NurbsCage    * pcage   , internal::gsXmlTree & data  );

  /// Extracts (parts of) a Brep from OpenNurbs
  bool readON_Brep        ( const ON_Brep         * pbrep   , internal::gsXmlTree & data  );

  /// Extracts a morph control object from OpenNurbs
  bool readON_MorphControl( const ON_MorphControl * pbrep   , internal::gsXmlTree & data  );
  
}

}
