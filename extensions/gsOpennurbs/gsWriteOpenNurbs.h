/** @file gsWriteOpenNurbs.h

    @brief Provides declaration of functions that write 3DM file format.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsDebug.h>

#include <gsIO/gsXmlUtils.h>

class ON_Surface;
class ON_Curve;
class ON_NurbsSurface;
class ON_NurbsCurve;
class ON_NurbsCage;
class ON_Brep;
class ON_MorphControl;


namespace gismo {

namespace extensions {

/// \todo Complete writing to onurbs
/*
template <typename T>
void gsWriteOpenNurbs(const gsGeometry<T>& geom,
                      const std::string& fileName);
*/

/// Writes a planar domain to OpenNurbs file
template<class T>
bool writeON_PlanarDomain( const gsPlanarDomain<T> & pd);


// Writes a NURBS curve to OpenNurbs file
//template<class T>
//bool writeON_NurbsCurve( const gsCurve<T> & curve)

}

}
