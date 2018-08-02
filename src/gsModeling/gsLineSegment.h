/** @file gsLineSegment.h

    @brief Provides declaration of Line class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Groiss, J. Vogl, A. Mantzaflaris
*/

#pragma once

namespace gismo
{

/// Represents a line segment in \a d dimensions
template<int d, class T>
class gsLineSegment
{
public:

    /// Constructs a line segment starting at \a point1 and ending at \a point2
    gsLineSegment(const gsPoint<d, T>& point1, const gsPoint<d, T>& point2)
    : m_start(point1), m_end(point2) { }

    /**
     * @brief Tells wheter line intersects segment
     *
     * This function returns true iff a line between origin and end
     * is intersected by this line.
     *
     * @param[in] origin
     * @param[in] end
     *
     * @return bool value indicating whether segment is intersected or not
     */
    bool intersectSegment(const gsPoint<d, T>& origin, const gsPoint<d, T>& end)
    {
        /*
          // Alternative implemantation
          gsLineSegment<2, T> segmentLine(origin, end);

          gsMatrix<T, 2, 2> matrix;
        matrix.col(0) = m_direction;
        matrix.col(1) = segmentLine.m_direction;
        if (math::abs(matrix.determinant()) <= 1e-8) {
            gsWarn << "Lines are parallel or ident.\n";
            return false;
        }
        gsVector<T, 2> parameters = matrix.partialPivLu().solve(m_start - segmentLine.m_start);
        double iparam = parameters(1);
        return ( -1e-8 <= iparam && iparam <= 1 + 1e-8 );
        */

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

        const gsVector<T, d> dir = direction();
        gsVector<T, d> p = origin - m_start;
        T d1 = dir[0] * p[1] - dir[1] * p[0];
        p = (end - m_start);
        T d2 = dir[0] * p[1] - dir[1] * p[0];
        size_t i1 = d1 > 0 ? 2 : (d1 < 0 ? 1 : 0);
        size_t i2 = d2 > 0 ? 2 : (d2 < 0 ? 1 : 0);
        return  (i1 ^ i2)!=0;
    }

    inline const gsVector<T, d> direction() const { return m_end - m_start; }

    T length() const { return (m_end - m_start).norm(); }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
    gsPoint<d, T> m_start;
    gsPoint<d, T> m_end;

}; // class gsLineSegment


} // namespace gismo
