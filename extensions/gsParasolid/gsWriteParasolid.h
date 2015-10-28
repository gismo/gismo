/** @file gsWriteParasolid.h

    @brief Provides declaration of gsWriteParasolid functions.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/


#pragma once

#include <string>

#include <gsCore/gsForwardDeclarations.h>

#include <gsParasolid/gsPKSession.h>

typedef int PK_GEOM_t;
typedef int PK_BSURF_t;
typedef int PK_BCURVE_t;
typedef int PK_BODY_t;
typedef int PK_ASSEMBLY_t;
struct PK_UVBOX_s;

namespace gismo {

namespace extensions {
    
    /// Writes a gsSurface to a parasolid file
    /// \param gssurf a surface
    /// \param fname filename (without extension)
    template<class T>
    bool gsWriteParasolid( const gsGeometry<T> & ggeo, std::string const & filename );

    template<class T>
    bool gsWriteParasolid( const gsMultiPatch<T> & gssurfs, std::string const & filename );

    template<class T>
    bool gsWriteParasolid( const gsMesh<T>& mesh, std::string const & filename );

    template<class T>
    bool gsWriteParasolid( const gsTHBSpline<2, T>& thb, std::string const & filename );

    /// Translates a gsTensorBSpline to a PK_BSURF_t
    /// \param[in] bsp B-spline surface
    /// \param[out] bsurf Parasolid spline surface
    template<class T, class KnotVectorType> 
    bool createPK_BSURF( const gsTensorBSpline< 2, T, KnotVectorType > & bsp, PK_BSURF_t & bsurf,
			 bool closed_u = false, bool closed_v = false );

    /// Translates a gsBSpline to a PK_BCURVE_t
    /// \param[in] curve B-Spline surve
    /// \param[out] bcurve Parasolid spline curve
    template<class T> 
    bool createPK_BCURVE( const gsBSpline<T>& curve, PK_BCURVE_t& bcurve );

    /// Translates a gsGeometry to a PK_GEOM_t
    /// \param[in] ggeo inpute G+SMO geometry
    /// \param[out] pgeo Parasolid geometric entity
    template<class T> 
    bool createPK_GEOM( const gsGeometry<T> & ggeo, PK_GEOM_t & pgeo );

    /// Translates a gsMesh to PK_BODY_t
    /// \param[in] mesh input G+Smo mesh
    /// \param[out] body Parasolid wire body
    template<class T> 
    bool exportMesh( const gsMesh<T>& mesh, PK_BODY_t& body );

    /// Translates a THB-Spline surface to PK_BODY_t
    /// \param[in] surface THB-Spline surface
    /// \param[out] body Parasolid body
    template<class T> 
    bool exportTHBsurface( const gsTHBSpline<2, T>& surface, PK_ASSEMBLY_t& body );


}//extensions

}//gismo

