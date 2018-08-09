/** @file gsReadBrep.cpp

    @brief Reading OpenCascade .brep files

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsDebug.h>

#include <gsIO/gsXml.h>

#define gsEndl std::endl

namespace opencascade { template<class O> class handle; }
class TopoDS_Shape;
class Geom_Surface;
class Geom_BSplineSurface;
class Geom2d_Curve;
class Geom2d_BSplineCurve;

namespace gismo {

namespace extensions {


/// Reads a brep file
bool gsReadBrep(const char * filename, internal::gsXmlTree & data);

/// Extracts a TopoDS_Shape
bool readTopoDS_Shape( const TopoDS_Shape & shape, internal::gsXmlTree & data  );

/// Extracts a surface
bool readGeom_Surface( const opencascade::handle<Geom_Surface> & S, internal::gsXmlTree & data  );

/// Extracts a B-spline surface
bool readGeom_BSplineSurface( const opencascade::handle<Geom_BSplineSurface> & S, internal::gsXmlTree & data  );

/// Extracts a 2d curve
bool readGeom2d_Curve( const opencascade::handle<Geom2d_Curve> & C, internal::gsXmlTree & data  );

/// Extracts a 2d B-spline curve
bool readGeom2d_BSplineCurve( const opencascade::handle<Geom2d_BSplineCurve> & bsp2d, internal::gsXmlTree & data  );

}

}
