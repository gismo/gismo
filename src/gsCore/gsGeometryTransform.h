/** @file gsGeometryTransform.h

    @brief Implements a Geometry together with an affine transform.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsAffineFunction.h>

namespace gismo {

/**
   @brief Representation of a transformed geometry

    \ingroup function
    \ingroup Core
 */
template <typename T>
class gsGeometryTransform : public gsFunction<T>
{
protected:
    const gsGeometry<T> * m_geo;
    gsAffineFunction<T> m_af;
public:
    /// Shared pointer for gsGeometryTransform
    typedef memory::shared_ptr< gsGeometryTransform > Ptr;

    /// Unique pointer for gsGeometryTransform
    typedef memory::unique_ptr< gsGeometryTransform > uPtr;

    GISMO_CLONE_FUNCTION(gsGeometryTransform)

    /**
     * @brief Affine maps are the composition of a linear map with a translation
     *        this constructor takes the two components explicitly
     * @param geo pointer to the original geometry map
     * @param mat the matrix corresponding to the linear map
     * @param trans the vector corresponding to the translation
     */
    gsGeometryTransform(const gsGeometry<T> * geo,
                        const gsMatrix<T> mat, const gsVector<T> trans)
    : m_geo(geo), m_af(mat,trans) { }

    /**
     * @brief Affine maps are the composition of a linear map with a translation
     *        this constructor takes the two components explicitly
     * @param mat the matrix corresponding to the linear map
     * @param trans the vector corresponding to the translation
     */
    static uPtr make(const gsMatrix<T> mat, const gsVector<T> trans)
    { return uPtr(new gsGeometryTransform(mat, trans)); }

    virtual short_t domainDim() const { return m_geo->domainDim(); }
    virtual short_t targetDim() const { return m_geo->targetDim(); }

    gsMatrix<T> support() const { return m_geo->support(); }
    
    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    
    virtual void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    
    virtual void deriv2_into( const gsMatrix<T>& u, gsMatrix<T>& result ) const
    { GISMO_NO_IMPLEMENTATION }
};

template<class T> void
gsGeometryTransform<T>::eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_geo->eval_into(u, result);
    result = ( m_af.matrix() * result ).colwise() + m_af.translation();
}

template<class T> void
gsGeometryTransform<T>::deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    m_geo->deriv_into(u, result);
    auto Jt = result.reshaped( u.rows(), result.size()/u.rows() );
    gsMatrix<T> tmp;
    index_t n = m_geo->targetDim();
    for(index_t i = 0; i!=u.cols(); ++i)
    {
        tmp = Jt.middleCols(n*i,n);
        Jt.middleCols(n*i,n).noalias() = tmp * m_af.matrix().transpose();
    }
}

} // namespace gismo

//#ifndef GISMO_BUILD_LIB
//#include GISMO_HPP_HEADER(gsGeometryTransform.hpp)
//#endif
