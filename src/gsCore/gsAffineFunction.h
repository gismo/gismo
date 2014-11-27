/** @file gsAffineFunction.h

    @brief


    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
    Created on: 2014-11-27
*/
#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsFunction.h>
#pragma once


namespace gismo {

/**
 * representation of an affine function
 */
template <typename T>
class gsAffineFunction : public gsFunction<T>
{
protected:
    gsMatrix<T> m_mat;
    gsVector<T> m_trans;
public:
    gsAffineFunction(const gsMatrix<T> mat, const gsVector<T> trans);
    virtual int domainDim() const;
    virtual int targetDim() const;
    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    virtual void eval_component_into(const gsMatrix<T>& u,
                                     const index_t comp,
                                     gsMatrix<T>& result) const;
    virtual void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    virtual void deriv2_into( const gsMatrix<T>& u, gsMatrix<T>& result ) const;
};


} // namespace gismo

#ifndef GISMO_HEADERS_ONLY
#include GISMO_HPP_HEADER(gsAffineFunction.hpp)
#endif
