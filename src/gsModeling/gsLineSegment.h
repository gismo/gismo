/** @file gsLine.h

    @brief Provides declaration of Line class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Groiss, J. Vogl, A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsDebug.h>

namespace gismo
{
template<int dim, class T>
class gsLineSegment
{
public:
    gsLineSegment() {}
    gsLineSegment(const gsPoint<dim, T>& point1, const gsPoint<dim, T>& point2) : m_point(point1),
        m_direction((point2 - point1)) { }

    /**
     * @brief Tells wheter line intersects segment
     *
     * This function returns if a line between origin and end will be intersected by the line where this function
     * is called. Or with other words, if origin and end are on different sides of line.
     * @param[in] origin
     * @param[in] end
     *
     * @return bool value indicating wheter segment is intersected or not
     */
    bool intersectSegment(const gsPoint<dim, T>& origin, const gsPoint<dim, T>& end)
    {
        /*gsLineSegment<2, real_t> segmentLine(origin, end);

        gsMatrix<T, 2, 2> matrix;
        matrix.col(0) = m_direction;
        matrix.col(1) = segmentLine.m_direction;
        if (math::abs(matrix.determinant()) <= 1e-8) {
            gsWarn << "Lines are parallel or ident.\n";
            return false;
        }
        gsVector<T, 2> parameters = matrix.partialPivLu().solve(m_point - segmentLine.m_point);
        double iparam = parameters(1);
        return ( -1e-8 <= iparam && iparam <= 1 + 1e-8 );*/

        /*
         * p1: -1 0 1
         * p2: -1 0 1
         *
         * -1 -1 => false
         * -1  0 => true
         * -1  1 => true
         *  0 -1 => true
         *  0  0 => false
         *  0  1 => true
         *  1 -1 => true
         *  1  0 => true
         *  1  1 => false
         */

        gsVector<T, 2> p = origin - m_point;
        real_t d1 = m_direction[0] * p[1] - m_direction[1] * p[0];
        p = (end - m_point);
        real_t d2 = m_direction[0] * p[1] - m_direction[1] * p[0];
        size_t i1 = d1 > 0 ? 2 : d1 < 0 ? 1 : 0;
        size_t i2 = d2 > 0 ? 2 : d2 < 0 ? 1 : 0;
        return  (i1 ^ i2)!=0;
    }

    const gsVector<T, dim> & direction() const { return m_direction; }

    T length() const { return m_direction.norm(); }

private:
    gsPoint<dim, T> m_point;
    gsVector<T, dim> m_direction;
}; // class gsLineSegment

typedef gsLineSegment<2, real_t> gsLineSegment2D;
//typedef gsLineSegment<3, real_t> gsSegment3D;

} // namespace gismo
