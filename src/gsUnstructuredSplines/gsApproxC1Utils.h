/** @file gsApproxC1Utils.h

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/


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

template <class T>
class gsTraceBasis : public gismo::gsFunction<T>
{

protected:
    gsGeometry<T> & _geo;

    gsBasis<T> & m_basis_plus;
    gsBasis<T> & m_basis_geo;
    gsBSpline<T> & _m_basis_beta;

    mutable gsMapData<T> _tmp;

    bool m_isboundary;
    const index_t m_bfID, m_uv;


public:
    /// Shared pointer for gsTraceBasis
    typedef memory::shared_ptr< gsTraceBasis > Ptr;

    /// Unique pointer for gsTraceBasis
    typedef memory::unique_ptr< gsTraceBasis > uPtr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    gsTraceBasis(gsGeometry<T> & geo,
                 gsBasis<T> & basis_plus,
                 gsBasis<T> & basis_geo,
                 gsBSpline<T> & basis_beta,
                 bool isboundary,
                 const index_t bfID,
                 const index_t uv) :
            _geo(geo), m_basis_plus(basis_plus), m_basis_geo(basis_geo), _m_basis_beta(basis_beta),
            m_isboundary(isboundary), m_bfID(bfID), m_uv(uv), _traceBasis_piece(nullptr)
    {
        //_tmp.flags = NEED_JACOBIAN;
    }

    ~gsTraceBasis() { delete _traceBasis_piece; }

GISMO_CLONE_FUNCTION(gsTraceBasis)

    short_t domainDim() const {return 2;}

    short_t targetDim() const {return 1;}

    mutable gsTraceBasis<T> * _traceBasis_piece; // why do we need this?

    const gsFunction<T> & piece(const index_t k) const
    {
        delete _traceBasis_piece;
        _traceBasis_piece = new gsTraceBasis(_geo, m_basis_plus, m_basis_geo, _m_basis_beta,
                                             m_isboundary, m_bfID, m_uv);
        return *_traceBasis_piece;
    }

    // Input is parametric coordinates of 2-D \a mp
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize( targetDim() , u.cols() );

        // tau/p
        gsBSplineBasis<T> bsp_temp = dynamic_cast<gsBSplineBasis<> & >(m_basis_geo);

        real_t p = bsp_temp.degree();
        real_t tau_1 = bsp_temp.knots().at(p + 1); // p + 2

        gsMatrix<T> beta, N_0, N_1, N_i_plus, der_N_i_plus;

        if (!m_isboundary)
            _m_basis_beta.eval_into(u.row(m_uv),beta); // 1-dir == PatchID
        else
            beta.setZero(1, u.cols());

        m_basis_geo.evalSingle_into(0,u.row(1-m_uv),N_0); // u
        m_basis_geo.evalSingle_into(1,u.row(1-m_uv),N_1); // u

        m_basis_plus.evalSingle_into(m_bfID,u.row(m_uv),N_i_plus); // v
        m_basis_plus.derivSingle_into(m_bfID,u.row(m_uv),der_N_i_plus);

        gsMatrix<T> temp = beta.cwiseProduct(der_N_i_plus);
        result = N_i_plus.cwiseProduct(N_0 + N_1) - temp.cwiseProduct(N_1) * tau_1 / p;
    }

};


template <class T>
class gsNormalDerivBasis : public gismo::gsFunction<T>
{

protected:
    gsGeometry<T> & _geo;

    gsBasis<T> & m_basis_minus;
    gsBasis<T> & m_basis_geo;
    gsBSpline<T> & m_basis_alpha;

    mutable gsMapData<T> _tmp;

    bool m_isboundary;
    const index_t m_bfID, m_uv;


public:
    /// Shared pointer for gsNormalDerivBasis
    typedef memory::shared_ptr< gsNormalDerivBasis > Ptr;

    /// Unique pointer for gsNormalDerivBasis
    typedef memory::unique_ptr< gsNormalDerivBasis > uPtr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    gsNormalDerivBasis(gsGeometry<T> & geo,
                 gsBasis<T> & basis_minus,
                 gsBasis<T> & basis_geo,
                 gsBSpline<T> & basis_alpha,
                 bool isboundary,
                 const index_t bfID,
                 const index_t uv) :
            _geo(geo), m_basis_minus(basis_minus), m_basis_geo(basis_geo), m_basis_alpha(basis_alpha),
            m_isboundary(isboundary), m_bfID(bfID), m_uv(uv), _normalDerivBasis_piece(nullptr)
    {
        //_tmp.flags = NEED_JACOBIAN;
    }

    ~gsNormalDerivBasis() { delete _normalDerivBasis_piece; }

GISMO_CLONE_FUNCTION(gsNormalDerivBasis)

    short_t domainDim() const {return 2;}

    short_t targetDim() const {return 1;}

    mutable gsNormalDerivBasis<T> * _normalDerivBasis_piece; // why do we need this?

    const gsFunction<T> & piece(const index_t k) const
    {
        delete _normalDerivBasis_piece;
        _normalDerivBasis_piece = new gsNormalDerivBasis(_geo, m_basis_minus, m_basis_geo, m_basis_alpha,
                                             m_isboundary, m_bfID, m_uv);
        return *_normalDerivBasis_piece;
    }

    // Input is parametric coordinates of 2-D \a mp
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize( targetDim() , u.cols() );

        // tau/p
        gsBSplineBasis<T> bsp_temp = dynamic_cast<gsBSplineBasis<> & >(m_basis_geo);

        real_t p = bsp_temp.degree();
        real_t tau_1 = bsp_temp.knots().at(p + 1); // p + 2

        gsMatrix<T> alpha, N_1, N_j_minus;

        if (!m_isboundary)
            m_basis_alpha.eval_into(u.row(m_uv),alpha); // 1-dir == PatchID
        else
            alpha.setOnes(1, u.cols());

        m_basis_geo.evalSingle_into(1,u.row(1-m_uv),N_1); // u

        m_basis_minus.evalSingle_into(m_bfID,u.row(m_uv),N_j_minus); // v

        if (!m_isboundary)
            result = (m_uv == 0 ? -1 : 1) * alpha.cwiseProduct(N_j_minus.cwiseProduct(N_1)) * tau_1 / p;
        else
            result = (m_uv == 0 ? -1 : 1) * alpha.cwiseProduct(N_j_minus.cwiseProduct(N_1));
    }

};

}