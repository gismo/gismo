/** @file gsClosestPoint.h

    @brief Provides declarations of functions for determining the
    closest point and parameter on a geometry to a given point using
    Parasolid functionality.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): D. Mokris
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsParasolid/gsPKSession.h>

// TODO: Why can't I include gsFrustrum.h already here?

struct PK_VECTOR_s;
typedef struct PK_VECTOR_s PK_VECTOR_t;

namespace gismo
{
    namespace extensions
    {
        /// Converts a G+Smo vector to a Parasolid vector.
        template <class T>
        PK_VECTOR_t gsPK_VECTOR(const gsVector<T, 3>& gsVector);

        /**
           Finds the point on \a gsBSurf that is closest to \a point
           and writes its parameter values into \a result.
           Returns Parasolid error status.
        */
        template <class T>
        bool gsClosestParam(const gsTensorBSpline<2, T>& gsBSurf,
                            const gsVector<T, 3>& gsPoint,
                            gsVector<T, 2>& gsResult);

        /**
           For each three-dimensional point represented as a row of \a
           gsPoints finds the closest point on \a gsBSurf and writes
           its parameters as a column of \a gsResults. (This might
           sound strange but it corresponds to the ordering in
           gsFitting.  Returns Parasolid error status.
        */
        template <class T>
        bool gsClosestParam(const gsTensorBSpline<2, T>& gsBSurf,
                            const gsMatrix<T>& gsPoints,
                            gsMatrix<T>& gsResults);

        /**
           Finds the point on \a gsBSurf that is closest to \a point
           and writes it into \a result.
           Returns Parasolid error status.

           TODO: Not implemented, yet.
        */
        template <class T>
        bool gsClosestPoint(const gsTensorBSpline<2, T>& gsBSurf,
                            const gsVector<T, 3>& gsPoint,
                            gsVector<T, 3>& gsResult)
        {
            gsWarn << "This function is not implemented, yet." << std::endl;
            return false;
        }

    } // namespace extensions

} // namespace gismo
    
