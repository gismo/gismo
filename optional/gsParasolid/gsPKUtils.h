/** @file gsPKUtils.h

    @brief Provides headers of utilities functions used by other parts of gsParasolid.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris, D. Mokris
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>

#include <gsParasolid/gsPKSession.h>

typedef int PK_BSURF_t;
typedef int PK_BCURVE_t;
typedef int PK_GEOM_t;

namespace gismo {

	namespace extensions {

		// returns true la all multiplicities in muls are less (<) than deg + 1
        // parasolid restriction
		bool validMultiplicities(const std::vector<index_t>& mult,
								 const int deg);

		/// Translates a gsTensorBSpline to a PK_BSURF_t
		/// \param[in] bsp B-spline surface
		/// \param[out] bsurf Parasolid spline surface
		template<class T>
			bool createPK_BSURF( const gsTensorBSpline< 2,T> & bsp, PK_BSURF_t & bsurf,
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

	} // namespace extensions

} // namespace gismo
