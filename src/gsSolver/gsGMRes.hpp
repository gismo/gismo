/** @file gsGMRes.hpp

    @brief Preconditioned iterative solver using the generalized minimal residual method.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

// TODO: Fix matrices sizes such that we don't resize on every iteration! (default can be 100 + 100 +...)

namespace gismo
{

template<class T>
bool gsGMRes<T>::initIteration( const typename gsGMRes<T>::VectorType& rhs,
                                typename gsGMRes<T>::VectorType& x )
{
    if (Base::initIteration(rhs,x))
        return true;

    m_mat->apply(x,tmp);
    tmp = rhs - tmp;
    m_precond->apply(tmp, residual);
    beta = residual.norm(); // This is  ||r||

    m_error = beta/m_rhs_norm;
    if(m_error < m_tol)
        return true;

    v.push_back(residual/beta);
    g.setZero(2,1);
    g(0,0) = beta;
    Omega = gsMatrix<T>::Identity(2, 2);
    Omega_prev = gsMatrix<T>::Identity(2, 2);

    return false;
}

template<class T>
void gsGMRes<T>::finalizeIteration( typename gsGMRes<T>::VectorType& x )
{
    //Remove last row of H and g
    H.resize(m_num_iter,m_num_iter);
    H = H_prev.block(0,0,m_num_iter,m_num_iter);
    g_tmp.resize(m_num_iter,1);
    g_tmp = g.block(0,0,m_num_iter,1);

    //Solve H*y = g;
    solveUpperTriangular(H, g_tmp);

    //Create the matrix from the column matrix in v.
    gsMatrix<T> V(m_mat->rows(),m_num_iter);
    for (index_t k = 0; k< m_num_iter; ++k)
    {
        V.col(k) = v[k];
    }
    //Update solution
    x += V*y;

    // cleanup temporaries
    tmp.clear();
    g.clear();
    g_tmp.clear();
    h_tmp.clear();
    y.clear();
    w.clear();
    residual.clear();
    H_prev.clear();
    H.clear();
    Omega.clear();
    Omega_prev.clear();
    Omega_tmp.clear();
    Omega_prev_tmp.clear();
    v.clear();
}

template<class T>
bool gsGMRes<T>::step( typename gsGMRes<T>::VectorType& )
{
    // The iterate x is never updated! Use finalizeIteration to obtain x.
    const index_t k = m_num_iter-1;
    H.setZero(k+2,k+1);
    h_tmp.setZero(k+2,1);

    if (k != 0)
    {
        H.block(0,0,k+1,k) = H_prev;
    }

    Omega = gsMatrix<T>::Identity(k+2, k+2);
    m_mat->apply(v[k],tmp);
    m_precond->apply(tmp, w);

    for (index_t i = 0; i< k+1; ++i)
    {
        h_tmp(i,0) = (w.transpose()*v[i]).value(); //Typo h_l,k
        w = w - h_tmp(i,0)*v[i];
    }
    h_tmp(k+1,0) = w.norm();

  //  if (math::abs(h_tmp(k+1,0)) < 1e-16) //If exact solution
  //      return true;

    v.push_back(w/h_tmp(k+1,0));

    h_tmp = Omega_prev*h_tmp;
    H.block(0,k,k+2,1) = h_tmp;

    //Find coef in rotation matrix
    T sk = H(k+1,k)/(math::sqrt(H(k,k)*H(k,k) + H(k+1,k)*H(k+1,k)));
    T ck = H(k,  k)/(math::sqrt(H(k,k)*H(k,k) + H(k+1,k)*H(k+1,k)));
    Omega(k,k)   = ck; Omega(k,k+1)   = sk;
    Omega(k+1,k) =-sk; Omega(k+1,k+1) = ck;

    //Rotate H and g
    H = Omega*H;
    H_prev = H;
    g_tmp.setZero(k+2,1);
    if (k != 0)
        g_tmp.block(0,0,k+1,1) = g;
    else
        g_tmp = g;
    g.noalias() = Omega * g_tmp;

    T residualNorm2 = g(k+1,0)*g(k+1,0);
    m_error = math::sqrt(residualNorm2) / m_rhs_norm;
    if(m_error < m_tol)
        return true;

    //Resize rotation product
    Omega_prev_tmp.setZero(k+3,k+3);
    Omega_prev_tmp.block(0,0,k+2,k+2) = Omega_prev;
    Omega_prev_tmp(k+2,k+2) = 1;

    Omega_tmp.setZero(k+3,k+3);
    Omega_tmp.block(0,0,k+2,k+2) = Omega;
    Omega_tmp(k+2,k+2) = 1;

    Omega_prev = Omega_tmp*Omega_prev_tmp;
    return false;
}

}
