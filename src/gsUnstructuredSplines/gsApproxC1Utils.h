/** @file gsApproxC1Utils.h

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/
#pragma once



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
        //delete _alpha_piece;
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
        //delete _beta_piece;
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
            D0  = ev.col(m_uv);
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
        //delete _traceBasis_piece;
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
        //delete _normalDerivBasis_piece;
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


template <class T>
class gsVertexBasis : public gismo::gsFunction<T>
{

protected:
    const gsGeometry<T>    & m_geo;
    std::vector<gsBSplineBasis<T>>       & m_basis_plus;
    std::vector<gsBSplineBasis<T>>       & m_basis_minus;
    std::vector<gsBSplineBasis<T>>       & m_basis_geo;
    std::vector<gsBSpline<T>>            & m_alpha;
    std::vector<gsBSpline<T>>            & m_beta;

    const real_t & m_sigma;
    const std::vector<bool> & m_kindOfEdge;

    const index_t m_bfID;

    mutable gsMapData<T> _tmp;

public:
    /// Shared pointer for gsVertexBasis
    typedef memory::shared_ptr< gsVertexBasis > Ptr;

    /// Unique pointer for gsVertexBasis
    typedef memory::unique_ptr< gsVertexBasis > uPtr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    gsVertexBasis(const gsGeometry<T>    & geo,
                  std::vector<gsBSplineBasis<T>>       & basis_plus,
                  std::vector<gsBSplineBasis<T>>       & basis_minus,
                  std::vector<gsBSplineBasis<T>>       & basis_geo,
                  std::vector<gsBSpline<T>>& alpha,
                  std::vector<gsBSpline<T>>& beta,
                  const real_t & sigma,
                  const std::vector<bool> & kindOfEdge,
                  const index_t bfID
            ) : m_geo(geo), m_basis_plus(basis_plus), m_basis_minus(basis_minus), m_basis_geo(basis_geo),
            m_alpha(alpha), m_beta(beta), m_sigma(sigma), m_kindOfEdge(kindOfEdge), m_bfID(bfID),
            _vertexBasis_piece(nullptr)
    {
        //_tmp.flags = NEED_JACOBIAN;
    }

    ~gsVertexBasis() { delete _vertexBasis_piece; }

    GISMO_CLONE_FUNCTION(gsVertexBasis)

    short_t domainDim() const {return 2;}

    short_t targetDim() const {return 1;}

    mutable gsVertexBasis<T> * _vertexBasis_piece; // why do we need this?

    const gsFunction<T> & piece(const index_t k) const
    {
        //delete _vertexBasis_piece;
        _vertexBasis_piece = new gsVertexBasis(m_geo, m_basis_plus, m_basis_minus, m_basis_geo, m_alpha,
                                               m_beta, m_sigma, m_kindOfEdge, m_bfID);
        return *_vertexBasis_piece;
    }

    // Input is parametric coordinates of 2-D \a mp
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.resize( targetDim() , u.cols() );

        result.setZero();


        // Computing the basis functions at the vertex
        gsMatrix<> Phi(6,6);
        Phi.setIdentity();

        Phi.row(1) *= m_sigma;
        Phi.row(2) *= m_sigma;
        Phi.row(3) *= m_sigma * m_sigma;
        Phi.row(4) *= m_sigma * m_sigma;
        Phi.row(5) *= m_sigma * m_sigma;

        // Computing c, c+ and c-
        // Point zero
        gsMatrix<> zero;
        zero.setZero(2,1);

        std::vector<gsMatrix<>> c_0, c_1;
        std::vector<gsMatrix < >> c_0_plus, c_1_plus, c_2_plus;
        std::vector<gsMatrix < >> c_0_plus_deriv, c_1_plus_deriv, c_2_plus_deriv;
        std::vector<gsMatrix < >> c_0_minus, c_1_minus;
        for (index_t i = 0; i < 2; i++) // i == 0 == u , i == 1 == v
        {
            gsMatrix<> b_0, b_1;
            gsMatrix<> b_0_plus, b_1_plus, b_2_plus;
            gsMatrix<> b_0_plus_deriv, b_1_plus_deriv, b_2_plus_deriv;
            gsMatrix<> b_0_minus, b_1_minus;

            gsBSplineBasis<T> bsp_temp = dynamic_cast<gsBSplineBasis<> & >(m_basis_geo[1-i]);
            real_t p = bsp_temp.degree();
            real_t h_geo = bsp_temp.knots().at(p + 1);

            m_basis_geo[1-i].evalSingle_into(0, u.row(i),b_0); // first
            m_basis_geo[1-i].evalSingle_into(1, u.row(i),b_1); // second

            m_basis_plus[i].evalSingle_into(0, u.row(i),b_0_plus);
            m_basis_plus[i].evalSingle_into(1, u.row(i),b_1_plus);
            m_basis_plus[i].evalSingle_into(2, u.row(i),b_2_plus);

            m_basis_plus[i].derivSingle_into(0, u.row(i),b_0_plus_deriv);
            m_basis_plus[i].derivSingle_into(1, u.row(i),b_1_plus_deriv);
            m_basis_plus[i].derivSingle_into(2, u.row(i),b_2_plus_deriv);

            m_basis_minus[i].evalSingle_into(0, u.row(i),b_0_minus);
            m_basis_minus[i].evalSingle_into(1, u.row(i),b_1_minus);

            c_0.push_back(b_0 + b_1);
            c_1.push_back((h_geo / p) * b_1);

            c_0_minus.push_back(b_0_minus + b_1_minus);
            c_1_minus.push_back(h_geo/ (p-1) * b_1_minus);

            gsMatrix<> der_b_1_plus_0, der2_b_1_plus_0, der2_b_2_plus_0;
            m_basis_plus[i].derivSingle_into(1, zero.row(i), der_b_1_plus_0);
            m_basis_plus[i].deriv2Single_into(1, zero.row(i), der2_b_1_plus_0);
            m_basis_plus[i].deriv2Single_into(2, zero.row(i), der2_b_2_plus_0);

            real_t factor_c_1_plus = 1/der_b_1_plus_0(0,0);
            real_t factor2_c_1_plus = -der2_b_1_plus_0(0,0)/(der_b_1_plus_0(0,0)*der2_b_2_plus_0(0,0));
            real_t factor_c_2_plus = 1/der2_b_2_plus_0(0,0);

            c_0_plus.push_back(b_0_plus + b_1_plus + b_2_plus);
            c_1_plus.push_back(factor_c_1_plus * b_1_plus + factor2_c_1_plus * b_2_plus);
            c_2_plus.push_back(factor_c_2_plus * b_2_plus );

            c_0_plus_deriv.push_back(b_0_plus_deriv + b_1_plus_deriv + b_2_plus_deriv);
            c_1_plus_deriv.push_back(factor_c_1_plus * b_1_plus_deriv + factor2_c_1_plus * b_2_plus_deriv);
            c_2_plus_deriv.push_back(factor_c_2_plus * b_2_plus_deriv);
        }

        std::vector<gsMatrix<>> alpha, beta, alpha_0, beta_0, alpha_deriv, beta_deriv;

        gsMatrix < T > temp_mat;
        if (m_kindOfEdge[0])
        {
            m_alpha[0].eval_into(u.row(0),temp_mat); // 1-dir == PatchID
            alpha.push_back(temp_mat); // u

            m_alpha[0].eval_into(zero.row(0),temp_mat); // 1-dir == PatchID
            alpha_0.push_back(temp_mat); // u

            m_alpha[0].deriv_into(zero.row(0),temp_mat); // 1-dir == PatchID
            alpha_deriv.push_back(temp_mat); // u

            m_beta[0].eval_into(u.row(0),temp_mat); // 1-dir == PatchID
            beta.push_back(temp_mat); // u

            m_beta[0].eval_into(zero.row(0),temp_mat); // 1-dir == PatchID
            beta_0.push_back(temp_mat); // u

            m_beta[0].deriv_into(zero.row(0),temp_mat); // 1-dir == PatchID
            beta_deriv.push_back(temp_mat); // u
        }
        else
        {
            temp_mat.setOnes(1, u.cols());
            alpha.push_back(temp_mat); // u

            temp_mat.setOnes(1, zero.cols());
            alpha_0.push_back(temp_mat); // u

            temp_mat.setZero(1, zero.cols());
            alpha_deriv.push_back(temp_mat); // u

            temp_mat.setZero(1, u.cols());
            beta.push_back(temp_mat); // u

            temp_mat.setZero(1, zero.cols());
            beta_0.push_back(temp_mat); // u

            temp_mat.setZero(1, zero.cols());
            beta_deriv.push_back(temp_mat); // u
        }



        if (m_kindOfEdge[1]) {
            m_alpha[1].eval_into(u.row(1), temp_mat); // 1-dir == PatchID
            alpha.push_back(temp_mat); // v

            m_alpha[1].eval_into(zero.row(0), temp_mat); // 1-dir == PatchID
            alpha_0.push_back(temp_mat); // v

            m_alpha[1].deriv_into(zero.row(0), temp_mat); // 1-dir == PatchID
            alpha_deriv.push_back(temp_mat); // v

            m_beta[1].eval_into(u.row(1), temp_mat); // 1-dir == PatchID
            beta.push_back(temp_mat); // v

            m_beta[1].eval_into(zero.row(0), temp_mat); // 1-dir == PatchID
            beta_0.push_back(temp_mat); // v

            m_beta[1].deriv_into(zero.row(0), temp_mat); // 1-dir == PatchID
            beta_deriv.push_back(temp_mat); // v
        }
        else
        {
            temp_mat.setOnes(1, u.cols());
            alpha.push_back(temp_mat); // u

            temp_mat.setOnes(1, zero.cols());
            alpha_0.push_back(temp_mat); // u

            temp_mat.setZero(1, zero.cols());
            alpha_deriv.push_back(temp_mat); // u

            temp_mat.setZero(1, u.cols());
            beta.push_back(temp_mat); // u

            temp_mat.setZero(1, zero.cols());
            beta_0.push_back(temp_mat); // u

            temp_mat.setZero(1, zero.cols());
            beta_deriv.push_back(temp_mat); // u
        }

        // Geo data:
        gsMatrix<> geo_jac = m_geo.jacobian(zero);
        gsMatrix<T> geo_der2 = m_geo.deriv2(zero);

        // Compute dd^^(i_k) and dd^^(i_k-1)
        gsMatrix<> dd_ik_plus, dd_ik_minus;
        gsMatrix<> dd_ik_minus_deriv, dd_ik_plus_deriv;
        dd_ik_minus = -1/(alpha_0[0](0,0)) * (geo_jac.col(1) +
                                              beta_0[0](0,0) * geo_jac.col(0));

        dd_ik_plus = 1/(alpha_0[1](0,0)) * (geo_jac.col(0) +
                                            beta_0[1](0,0) * geo_jac.col(1));

        gsMatrix<> geo_deriv2_12(2,1), geo_deriv2_11(2,1), geo_deriv2_22(2,1);
        geo_deriv2_12.row(0) = geo_der2.row(2);
        geo_deriv2_12.row(1) = geo_der2.row(5);
        geo_deriv2_11.row(0) = geo_der2.row(0);
        geo_deriv2_11.row(1) = geo_der2.row(3);
        geo_deriv2_22.row(0) = geo_der2.row(1);
        geo_deriv2_22.row(1) = geo_der2.row(4);
        gsMatrix<> alpha_squared_u = alpha_0[0]*alpha_0[0];
        gsMatrix<> alpha_squared_v = alpha_0[1]*alpha_0[1];

        dd_ik_minus_deriv = -1/(alpha_squared_u(0,0)) * // N^2
                            ((geo_deriv2_12 + (beta_deriv[0](0,0) * geo_jac.col(0) +
                                               beta_0[0](0,0) * geo_deriv2_11))*alpha_0[0](0,0) -
                             (geo_jac.col(1) + beta_0[0](0,0) * geo_jac.col(0)) *
                             alpha_deriv[0](0,0));

        dd_ik_plus_deriv = 1/(alpha_squared_v(0,0)) *
                           ((geo_deriv2_12 + (beta_deriv[1](0,0) * geo_jac.col(1) +
                                              beta_0[1](0,0) * geo_deriv2_22))*alpha_0[1](0,0) -
                            (geo_jac.col(0) + beta_0[1](0,0) * geo_jac.col(1)) *
                            alpha_deriv[1](0,0));

        // Comupute d_(0,0)^(i_k), d_(1,0)^(i_k), d_(0,1)^(i_k), d_(1,1)^(i_k) ; i_k == 2
        std::vector<gsMatrix<>> d_ik;
        d_ik.push_back(Phi.col(0));
        d_ik.push_back(Phi.block(0,1,6,2) * geo_jac.col(0) ); // deriv into u
        d_ik.push_back(Phi.block(0,1,6,2) * geo_jac.col(1) ); // deriv into v
        d_ik.push_back((geo_jac(0,0) * Phi.col(3) + geo_jac(1,0) * Phi.col(4))*geo_jac(0,1) +
                       (geo_jac(0,0) * Phi.col(4) + geo_jac(1,0) * Phi.col(5))*geo_jac(1,1) +
                       Phi.block(0,1,6,1) * geo_der2.row(2) +
                       Phi.block(0,2,6,1) * geo_der2.row(5)); // Hessian

        // Compute d_(*,*)^(il,ik)
        std::vector<gsMatrix<>> d_ilik_minus, d_ilik_plus;
        d_ilik_minus.push_back(Phi.col(0));
        d_ilik_minus.push_back(Phi.block(0,1,6,2) * geo_jac.col(0));
        d_ilik_minus.push_back((geo_jac(0,0) * Phi.col(3) + geo_jac(1,0) * Phi.col(4))*geo_jac(0,0) +
                               (geo_jac(0,0) * Phi.col(4) + geo_jac(1,0) * Phi.col(5))*geo_jac(1,0) +
                               Phi.block(0,1,6,1) * geo_der2.row(0) +
                               Phi.block(0,2,6,1) * geo_der2.row(3));
        d_ilik_minus.push_back(Phi.block(0,1,6,2) * dd_ik_minus);
        d_ilik_minus.push_back((geo_jac(0,0) * Phi.col(3) + geo_jac(1,0) * Phi.col(4))*dd_ik_minus(0,0) +
                               (geo_jac(0,0) * Phi.col(4) + geo_jac(1,0) * Phi.col(5))*dd_ik_minus(1,0) +
                               Phi.block(0,1,6,1) * dd_ik_minus_deriv.row(0) +
                               Phi.block(0,2,6,1) * dd_ik_minus_deriv.row(1));

        d_ilik_plus.push_back(Phi.col(0));
        d_ilik_plus.push_back(Phi.block(0,1,6,2) * geo_jac.col(1));
        d_ilik_plus.push_back((geo_jac(0,1) * Phi.col(3) + geo_jac(1,1) * Phi.col(4))*geo_jac(0,1) +
                              (geo_jac(0,1) * Phi.col(4) + geo_jac(1,1) * Phi.col(5))*geo_jac(1,1) +
                              Phi.block(0,1,6,1) * geo_der2.row(1) +
                              Phi.block(0,2,6,1) * geo_der2.row(4));
        d_ilik_plus.push_back(Phi.block(0,1,6,2) * dd_ik_plus);
        d_ilik_plus.push_back((geo_jac(0,1) * Phi.col(3) + geo_jac(1,1) * Phi.col(4))*dd_ik_plus(0,0) +
                              (geo_jac(0,1) * Phi.col(4) + geo_jac(1,1) * Phi.col(5))*dd_ik_plus(1,0) +
                              Phi.block(0,1,6,1) * dd_ik_plus_deriv.row(0) +
                              Phi.block(0,2,6,1) * dd_ik_plus_deriv.row(1));



        result = d_ilik_minus.at(0)(m_bfID,0) * (c_0_plus.at(0).cwiseProduct(c_0.at(1)) -
                                                   beta[0].cwiseProduct(c_0_plus_deriv.at(0).cwiseProduct(c_1.at(1)))) +
                        d_ilik_minus.at(1)(m_bfID,0) * (c_1_plus.at(0).cwiseProduct(c_0.at(1)) -
                                                   beta[0].cwiseProduct(c_1_plus_deriv.at(0).cwiseProduct(c_1.at(1)))) +
                        d_ilik_minus.at(2)(m_bfID,0) * (c_2_plus.at(0).cwiseProduct(c_0.at(1)) -
                                                   beta[0].cwiseProduct(c_2_plus_deriv.at(0).cwiseProduct(c_1.at(1)))) -
                        d_ilik_minus.at(3)(m_bfID,0) * alpha[0].cwiseProduct(c_0_minus.at(0).cwiseProduct(c_1.at(1))) -
                        d_ilik_minus.at(4)(m_bfID,0) * alpha[0].cwiseProduct(c_1_minus.at(0).cwiseProduct(c_1.at(1))); // f*_(ik-1,ik)

        //if (kindOfEdge[0])
        //rhsVals.at(i).setZero();

        //if (!kindOfEdge[1])
        result += d_ilik_plus.at(0)(m_bfID,0) * (c_0_plus.at(1).cwiseProduct(c_0.at(0)) -
                                                   beta[1].cwiseProduct(c_0_plus_deriv.at(1).cwiseProduct(c_1.at(0)))) +
                         d_ilik_plus.at(1)(m_bfID,0) * (c_1_plus.at(1).cwiseProduct(c_0.at(0)) -
                                                   beta[1].cwiseProduct(c_1_plus_deriv.at(1).cwiseProduct(c_1.at(0)))) +
                         d_ilik_plus.at(2)(m_bfID,0) * (c_2_plus.at(1).cwiseProduct(c_0.at(0)) -
                                                   beta[1].cwiseProduct(c_2_plus_deriv.at(1).cwiseProduct(c_1.at(0)))) +
                         d_ilik_plus.at(3)(m_bfID,0) * alpha[1].cwiseProduct(c_0_minus.at(1).cwiseProduct(c_1.at(0))) +
                         d_ilik_plus.at(4)(m_bfID,0) * alpha[1].cwiseProduct(c_1_minus.at(1).cwiseProduct(c_1.at(0))); // f*_(ik+1,ik)

        result -= d_ik.at(0)(m_bfID,0) * c_0.at(0).cwiseProduct(c_0.at(1)) + d_ik.at(2)(m_bfID,0) * c_0.at(0).cwiseProduct(c_1.at(1)) +
                         d_ik.at(1)(m_bfID,0) * c_1.at(0).cwiseProduct(c_0.at(1)) + d_ik.at(3)(m_bfID,0) * c_1.at(0).cwiseProduct(c_1.at(1)); // f*_(ik)


    }

};

}