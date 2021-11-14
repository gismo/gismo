/** @file gsApproxGluingData.hpp

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#include<gsUnstructuredSplines/gsApproxGluingData.h>

namespace gismo
{

// Input is parametric coordinates of 1-D \a mp
template <class T>
class gsAlpha : public gismo::gsFunction<T>
{

protected:
    gsGeometry<T> & _geo;
    mutable gsMapData<T> _tmp;
    index_t m_uv;


public:
    /// Shared pointer for gsAlpha
    typedef memory::shared_ptr< gsAlpha > Ptr;

    /// Unique pointer for gsAlpha
    typedef memory::unique_ptr< gsAlpha > uPtr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    gsAlpha(gsGeometry<T> & geo, index_t uv) :
            _geo(geo), m_uv(uv), _alpha_piece(nullptr)
    {
        _tmp.flags = NEED_JACOBIAN;
    }

    ~gsAlpha() { delete _alpha_piece; }

    GISMO_CLONE_FUNCTION(gsAlpha)

    short_t domainDim() const {return 1;}

    short_t targetDim() const {return 1;}

    mutable gsAlpha<T> * _alpha_piece; // why do we need this?

    const gsFunction<T> & piece(const index_t k) const
    {
        delete _alpha_piece;
        _alpha_piece = new gsAlpha(_geo, m_uv);
        return *_alpha_piece;
    }

    // Input is parametric coordinates of 1-D \a mp
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize( targetDim() , u.cols() );

        gsMatrix<T> uv, ev, D0;
        uv.setZero(2,u.cols());
        uv.row(m_uv) = u; // u

        T gamma = 1.0;

        for (index_t i = 0; i < uv.cols(); i++)
        {
            _geo.jacobian_into(uv.col(i), ev);
            uv(0, i) = gamma * ev.determinant();
        }
        result = uv.row(0);
    }

};

template <class T>
class gsBeta : public gismo::gsFunction<T>
{

protected:
    gsGeometry<T> & _geo;
    mutable gsMapData<T> _tmp;
    index_t m_uv;


public:
    /// Shared pointer for gsBeta
    typedef memory::shared_ptr< gsBeta > Ptr;

    /// Unique pointer for gsBeta
    typedef memory::unique_ptr< gsBeta > uPtr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    gsBeta(gsGeometry<T> & geo, index_t uv) :
            _geo(geo), m_uv(uv), _beta_piece(nullptr)
    {
        _tmp.flags = NEED_JACOBIAN;
    }

    ~gsBeta() { delete _beta_piece; }

    GISMO_CLONE_FUNCTION(gsBeta)

    short_t domainDim() const {return 1;}

    short_t targetDim() const {return 1;}

    mutable gsBeta<T> * _beta_piece; // why do we need this?

    const gsFunction<T> & piece(const index_t k) const
    {
        delete _beta_piece;
        _beta_piece = new gsBeta(_geo, m_uv);
        return *_beta_piece;
    }

    // Input is parametric coordinates of 1-D \a mp
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize( targetDim() , u.cols() );

        gsMatrix<T> uv, ev, D0;

        uv.setZero(2,u.cols());
        uv.row(m_uv) = u; // u

        T gamma = 1.0;

        for(index_t i = 0; i < uv.cols(); i++)
        {
            _geo.jacobian_into(uv.col(i),ev);
            D0 = ev.col(m_uv);
            real_t D1 = 1/ D0.norm();
            uv(0,i) = - gamma * D1 * D1 * ev.col(1).transpose() * ev.col(0);
        }
        result = uv.row(0);
    }

};

} // namespace gismo
