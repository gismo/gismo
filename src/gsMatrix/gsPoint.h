/** @file gsPoint.h

    @brief Provides declaration of Point class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Groiss, J. Vogl, A. Mantzaflaris
*/

#pragma once

// Assumes that Eigen library has been already included - like with gsVector.h

namespace gismo
{
/**
 * @brief A Point in T^dim, with an index number.
 * @tparam dim
 * @tparam T
 */
template<int dim, class T>
class gsPoint : public gsVector<T, dim>
{
public:
    typedef T Scalar_t;

    //gsPoint() : gsVector<T, dim>(), m_vertexIndex((size_t)0) {}

    //gsPoint(const gsPoint<dim, T>& point) : gsVector<T,dim>(), m_vertexIndex(point.m_vertexIndex) { *this << point[0], point[1]; }

    gsPoint(T x, T y, size_t index) : gsVector<T, 2>(), m_vertexIndex(index) { *this << x, y; }
    gsPoint(T x, T y, T z, size_t index) : gsVector<T, 3>(), m_vertexIndex(index) { *this << x, y, z; }

    inline int getVertexIndex() const { return m_vertexIndex; }

    void moveToPosition(const T x, const T y) {
        *this << x, y;
    }
    void moveToPosition(const T x, const T y, const T z) {
        *this << x, y, z;
    }

    void shift(const T x, const T y) {
        *this += gsVector<>::vec(x, y);
    }

    void shift(const T x, const T y, const T z) {
        *this += gsVector<>::vec(x, y, z);
    }

    static gsPoint<2,T> point(T x, T y, size_t index)
    {
        return gsPoint<2, T>(x, y, index);
    };

    static gsPoint<3,T> point(T x, T y, T z, size_t index)
    {
        return gsPoint<3, T>(x, y, z, index);
    };

private:
    size_t m_vertexIndex;
}; // class gsPoint

typedef gsPoint<2, real_t> gsPoint2D;
typedef gsPoint<3, real_t> gsPoint3D;

} // namespace gismo
