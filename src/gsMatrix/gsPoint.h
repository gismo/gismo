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
 * @brief A Point in T^d, with an index number.
 * @tparam d
 * @tparam T
 */
template<int d, class T>
class gsPoint : public gsVector<T, d>
{
public:
    typedef gsVector<T, d> Base;

    typedef gsPoint<d, T> Self;

    typedef typename Eigen::aligned_allocator<Self> aalloc;

    typedef T Scalar_t;

    gsPoint() : Base(), m_vertexIndex((size_t)0) {}

    gsPoint(T x, T y, size_t index) : Base(), m_vertexIndex(index) { *this << x, y; }

    inline int getVertexIndex() const { return m_vertexIndex; }

private:
    size_t m_vertexIndex;
}; // class gsPoint

} // namespace gismo
