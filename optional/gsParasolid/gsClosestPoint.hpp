/** @file gsClosestPoint.h

    @brief Provides definitions of functions for determining the
    closest point and parameter on a geometry to a given point using
    Parasolid functionality.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): D. Mokris
*/

#include <gsParasolid/gsFrustrum.h>
#include <gsParasolid/gsClosestPoint.h>
#include <gsParasolid/gsWriteParasolid.h>

namespace gismo
{
    namespace extensions
    {
        template <class T>
        PK_VECTOR_t gsPK_VECTOR(const gsVector<T, 3>& gsVector)
        {
            PK_VECTOR_t result;

            for(index_t i=0; i<3; i++)
                result.coord[i] = double(gsVector(i));

            return result;
        }

        template <class T>
        bool gsClosestParam(const gsTensorBSpline<2, T>& gsBSurf,
                            const gsVector<T, 3>& gsPoint,
                            gsVector<T, 2>& gsResult)
        {
            PK_BSURF_t               bsurf;
            createPK_BSURF<T>(gsBSurf, bsurf, false, false);
            // TODO: What shall we do about periodic surfaces?

            PK_VECTOR_t              point = gsPK_VECTOR(gsPoint);

            PK_GEOM_range_vector_o_t options;
            PK_GEOM_range_vector_o_m(options);
            options.opt_level =      PK_range_opt_accuracy_c;

            PK_range_result_t        range_result;
            PK_range_1_r_t           range;

            PK_ERROR_code_t err =    PK_GEOM_range_vector(bsurf, point, &options,
                                                          &range_result, &range);
            if(err)
                gsWarn << "err: " << err << std::endl;

            for(index_t i=0; i<2; i++)
                gsResult(i) = range.end.parameters[i];

            return err;
        }

        template <class T>
        bool gsClosestParam(const gsTensorBSpline<2, T>& gsBSurf,
                            const gsMatrix<T>& gsPoints,
                            gsMatrix<T>& gsResults)
        {
            GISMO_ASSERT(gsPoints.cols() == 3, "gsClosestParam is implemented for three-dimensional points only.");
            PK_BSURF_t               bsurf;
            createPK_BSURF<T>(gsBSurf, bsurf, false, false);

            int                      n_vectors = gsPoints.rows();
            PK_VECTOR_t              vectors[n_vectors];
            for(int i=0; i<n_vectors; i++)
            {
                gsVector<T, 3> point = gsPoints.row(i);
                vectors[i] = gsPK_VECTOR(point);
            }

            PK_GEOM_range_vector_many_o_t options;
            PK_GEOM_range_vector_many_o_m(options);

            // Keeping the opt_level at PK_range_opt_performace_c
            // (default) lead to slightly different results between
            // the two versions of this function.
            options.opt_level =      PK_range_opt_accuracy_c;

            PK_range_result_t        range_result[n_vectors];
            PK_range_1_r_t           ranges[n_vectors];

            PK_ERROR_code_t err =    PK_GEOM_range_vector_many(bsurf, n_vectors, vectors, &options,
                                                               range_result, ranges);
            if(err)
                gsWarn << "err: " << err << std::endl;

            // gsResults.col(i) = (u[i], v[i])
            gsResults.resize(2, n_vectors);
            for(int j=0; j<n_vectors; j++)
                for(int i=0; i<2; i++)
                    gsResults(i, j) = ranges[j].end.parameters[i];

            return err;
        }

    } // namespace extensions

} // namespace gismo
