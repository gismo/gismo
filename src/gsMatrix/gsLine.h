/** @file gsLine.h

    @brief Provides declaration of Line class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Groiss, J. Vogl, A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsDebug.h>

#define PRECISION 0.00000000001

namespace gismo
{
template<int dim, class T>
class gsLine : public gsPoint<dim, T>
{
public:
    gsLine() : gsPoint<dim, T>(), m_direction() {}
    gsLine(const gsPoint<dim, T>& point1, const gsPoint<dim, T>& point2) : gsPoint<dim, T>(point1),
        m_direction(point2 - point1) { }

    /**
     * @brief Calculates intersection parameter
     *
     * This function calculates the intersection parameter of two lines according to THIS.
     * The two lines are given like
     *  THIS: m_point + s*m_direction
     *  line2D: m_point + t*m_direction
     * and s is returned.
     *
     * In case the directions equal each other an error message is returned and 0 is returned.
     *
     * @param[in] line2D const line2D& - line
     * @return intersection parameter
     */
    double intersectionParameter(const gsLine<dim, T>& line) const
    {
        if(m_direction == line.m_direction)
        {
            gsWarn << "Lines are parallel or ident.\n";
            return 0;
        }
        else
        {
            gsMatrix<T,2,2> matrix;
            matrix(0,0) = line.m_direction[0];
            matrix(0,1) = m_direction[0];
            matrix(1,0) = line.m_direction[1];
            matrix(1,1) = m_direction[1];
            gsVector<T,2> rhs;
            rhs << (*this)(0) - line[0], (*this)(1) - line[1];
            gsVector<T,2> parameters = matrix.partialPivLu().solve(rhs);
            return -parameters(1);
        }
    }

    /**
     * @brief Tells wheter line intersects segment
     *
     * This function returns
     *  TRUE if -myConstants::PRECISION <= intersectionParameter < 1 + myConstants::PRECISION
     *  FALSE otherwise.
     *
     * @param[in] origin const Point2D& - point
     * @param[in] end const Point2D& - point
     *
     * @return bool value indicating wheter segment is intersected or not
     */
    bool intersectSegment(const gsPoint<dim, T>& origin, const gsPoint<dim, T>& end)
    {
        gsLine<2, real_t> segmentLine(origin, end);
        if (m_direction == segmentLine.m_direction)
            return false;
        else if ( -PRECISION <= segmentLine.intersectionParameter(*this) && segmentLine.intersectionParameter(*this) < 1 + PRECISION)
            return true;
        else
            return false;
    }

    bool sameDirection(const gsLine<dim, T>& line)
    {
        return m_direction == line.m_direction;
    }

private:
    gsVector<T, dim> m_direction;
    
}; // class gsLine

typedef gsLine<2, real_t> gsLine2D;
//typedef gsLine<3, real_t> gsLine3D;

} // namespace gismo
