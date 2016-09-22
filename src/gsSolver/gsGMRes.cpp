/** @file gsGMRes.cpp

    @brief Preconditioned iterative solver using the generalized minimal residual method.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/
#include <gsSolver/gsGMRes.h>

// TODO
// Fix matrices sizes such that we don't resize on every iteration! (default can be 100 + 100 +...)
//
namespace gismo
{

void gsGMRes::initIteration( const VectorType& rhs, const VectorType& x0, const gsLinearOperator<>& precond)
{
    GISMO_ASSERT(rhs.cols()== 1, "Implemented only for single columns right hand side matrix");
    m_rhs = rhs;
    xInit = x0;
    m_mat->apply(x0,tmp);
    tmp = m_rhs - tmp;
    precond.apply(tmp, residual);
    beta = residual.norm(); // This is  ||r||
    v.push_back(residual/beta);
    g.setZero(2,1);
    g(0,0) = beta;
    Omega = gsMatrix<real_t>::Identity(2, 2);
    Omega_prew = gsMatrix<real_t>::Identity(2, 2);

    rhsNorm2 = rhs.squaredNorm(); // This is ||r||²
    if (rhsNorm2 == 0)
        rhsNorm2 = 1.0;
    residualNorm2 = 0;
    threshold = m_tol*m_tol*rhsNorm2;
    m_num_iter = 0;
}

void gsGMRes::solve(const VectorType& rhs, VectorType& x, const gsLinearOperator<>& precond)
{
    initIteration(rhs, x, precond);

    while(m_num_iter < m_max_iters)
    {
        m_num_iter++;
        if (step(x, precond))
            break;
    }
    m_error = math::sqrt(residualNorm2 / rhsNorm2);

    //Post processing
    //Remove last row of H and g
    H.resize(m_num_iter,m_num_iter);
    H = H_prew.block(0,0,m_num_iter,m_num_iter);
    g_tmp.resize(m_num_iter,1);
    g_tmp = g.block(0,0,m_num_iter,1);

    //Solve H*y = g;
    solveUpperTriangular(H, g_tmp);

    //Create the matrix from the column matrix in v.
    gsMatrix<real_t> V(m_mat->rows(),m_num_iter);
    for (index_t k = 0; k< m_num_iter; ++k)
    {
        V.col(k) = v[k];
    }
    //Update solution
    x = xInit + V*y;
}

bool gsGMRes::step( VectorType& x, const gsLinearOperator<>& precond )
{
    GISMO_UNUSED(x);
    const index_t k = m_num_iter-1;
    H.setZero(k+2,k+1);
    h_tmp.setZero(k+2,1);

    if (k != 0)
    {
        H.block(0,0,k+1,k) = H_prew;
    }

    Omega = gsMatrix<real_t>::Identity(k+2, k+2);
    m_mat->apply(v[k],tmp);
    precond.apply(tmp, w);

    for (index_t i = 0; i< k+1; ++i)
    {
        h_tmp(i,0) = (w.transpose()*v[i]).value(); //Typo h_l,k
        w = w - h_tmp(i,0)*v[i];
    }
    h_tmp(k+1,0) = w.norm();

    if (math::abs(h_tmp(k+1,0)) < 1e-16) //If exact solution
        return true;

    v.push_back(w/h_tmp(k+1,0));

    h_tmp = Omega_prew*h_tmp;
    H.block(0,k,k+2,1) = h_tmp;

    //Find coef in rotation matrix
    real_t sk = H(k+1,k)/(math::sqrt(H(k,k)*H(k,k) + H(k+1,k)*H(k+1,k)));
    real_t ck = H(k,  k)/(math::sqrt(H(k,k)*H(k,k) + H(k+1,k)*H(k+1,k)));
    Omega(k,k)   = ck; Omega(k,k+1)   = sk;
    Omega(k+1,k) =-sk; Omega(k+1,k+1) = ck;

    //Rotate H and g
    H = Omega*H;
    H_prew = H;
    g_tmp.setZero(k+2,1);
    if (k != 0)
        g_tmp.block(0,0,k+1,1) = g;
    else
        g_tmp = g;
    g = Omega*g_tmp;

    residualNorm2 = g(k+1,0)*g(k+1,0);
    if(residualNorm2 < threshold)
        return true;

    //Resize rotation product
    Omega_prew_tmp.setZero(k+3,k+3);
    Omega_prew_tmp.block(0,0,k+2,k+2) = Omega_prew;
    Omega_prew_tmp(k+2,k+2) = 1;

    Omega_tmp.setZero(k+3,k+3);
    Omega_tmp.block(0,0,k+2,k+2) = Omega;
    Omega_tmp(k+2,k+2) = 1;

    Omega_prew = Omega_tmp*Omega_prew_tmp;
    return false;
}

}
