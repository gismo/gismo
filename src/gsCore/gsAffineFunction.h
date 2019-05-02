/** @file gsAffineFunction.h

    @brief Implements an affine function.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
    Created on: 2014-11-27
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsFunction.h>

namespace gismo {

/**
   @brief Representation of an affine function

    \ingroup function
    \ingroup Core
 */
template <typename T>
class gsAffineFunction : public gsFunction<T>
{
protected:
    gsMatrix<T> m_mat;
    gsVector<T> m_trans;
public:
    /// Shared pointer for gsAffineFunction
    typedef memory::shared_ptr< gsAffineFunction > Ptr;

    /// Unique pointer for gsAffineFunction
    typedef memory::unique_ptr< gsAffineFunction > uPtr;

    /**
     * @brief copy constructor
     */
    gsAffineFunction(const gsAffineFunction &other)
        : m_mat(other.m_mat), m_trans(other.m_trans)
    {}

    GISMO_CLONE_FUNCTION(gsAffineFunction)

    /**
     * @brief Affine maps are the composition of a linear map with a translation
     *        this constructor takes the two components explicitly
     * @param mat the matrix corresponding to the linear map
     * @param trans the vector corresponding to the translation
     */
    gsAffineFunction(const gsMatrix<T> mat, const gsVector<T> trans);

    /**
     * @brief Construct the affine map that maps \a box1 to \a box2 by
     *        mapping coordinate \a i of \a box1 to coordinate \a directions(i) of \a box2
     * @param directions specifies the permutation of the coordinate directions, i.e. the rotation
     * @param orientation specifies the sign of the destination
     * @param box1 description of the domain box: lower corner in the first column, upper corner in the second
     * @param box2 description of the target box
     */
    gsAffineFunction(const gsVector<index_t> &directions, const gsVector<bool> &orientation, const gsMatrix<T> &box1, const gsMatrix<T> &box2);

    /**
     * @brief Affine maps are the composition of a linear map with a translation
     *        this constructor takes the two components explicitly
     * @param mat the matrix corresponding to the linear map
     * @param trans the vector corresponding to the translation
     */
    static uPtr make(const gsMatrix<T> mat, const gsVector<T> trans)
    { return uPtr(new gsAffineFunction(mat, trans)); }

    /**
     * @brief Construct the affine map that maps \a box1 to \a box2 by
     *        mapping coordinate \a i of \a box1 to coordinate \a directions(i) of \a box2
     * @param directions specifies the permutation of the coordinate directions, i.e. the rotation
     * @param orientation specifies the sign of the destination
     * @param box1 description of the domain box: lower corner in the first column, upper corner in the second
     * @param box2 description of the target box
     */
    static uPtr make(const gsVector<index_t> &directions, const gsVector<bool> &orientation, const gsMatrix<T> &box1, const gsMatrix<T> &box2)
    { return uPtr(new gsAffineFunction(directions, orientation, box1, box2)); }

    virtual short_t domainDim() const;
    virtual short_t targetDim() const;
    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    virtual void eval_component_into(const gsMatrix<T>& u,
                                     const index_t comp,
                                     gsMatrix<T>& result) const;
    virtual void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    virtual void deriv2_into( const gsMatrix<T>& u, gsMatrix<T>& result ) const;
};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsAffineFunction.hpp)
#endif
