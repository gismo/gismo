/** @file gsMinResQLP.hpp

    @brief Preconditioned iterative solver using the minimal residual QLP method.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

namespace gismo
{

//A stable way of finding radius, cosinus and sinus (r,c,s) from a and b
//That is: r = sqrt(a*a + b*b), c = a/r, s= b/r
//Rewritten from PETSc minres
template<class T>
void MinResQLPSymOrtho(T const a, T const b, T &c, T &s, T &r)
{
    if (b == 0.0) 
    {
        if (a == 0.0) 
            c = 1.0;
        else 
            c = math::getSign(a); 
        s = 0.0;
        r = math::abs(a); 
    } 
    else if (a == 0.0) 
    {
        c = 0.0;
        s = math::getSign(b);
        r = math::abs(b);
    } 
    else if (math::abs(b) > math::abs(a)) 
    {
        T t = a / b;

        s = math::getSign(b) / math::sqrt(1.0 + t * t);
        c = s * t;
        r = b / s;
    } 
    else
    {
        T t = b / a;
        c = math::getSign(a) / math::sqrt(1.0 + t * t);
        s = c * t;
        r = a / c;
    }
}


template<class T>
bool gsMinResQLP<T>::initIteration( const typename gsMinResQLP<T>::VectorType& rhs, typename gsMinResQLP<T>::VectorType& x )
{
    m_rhs = rhs;
    Base::initIteration(rhs,x);
    index_t n = m_mat->cols();
    index_t m = 1; // = rhs.cols();

    resvec.setZero(m_max_iters+1,1);
    Aresvec.setZero(m_max_iters+1,1);

    //r2 = rhs; 
    m_mat->apply(x,r2); //erex
    xnorm = x.norm(); 
    Axnorm = r2.norm();
    r2 = rhs-r2;
    //y = x + beta y.
    //VecAYPX(Vec y, PetscScalar beta, Vec x)
    r3.setZero(n,m); 
    m_precond->apply(r2, r3);
    beta1 = r2.col(0).dot(r3.col(0));
    GISMO_ASSERT(beta1 > 0, "Preconditioner is indefinite");
    beta1 = math::sqrt(beta1);
    beta  = 0; 
    tau   = 0; taul = 0;  gmin = 0; 
    phi   = beta1;  betan = beta1;
    cs = -1; cr1 = -1; cr2 = -1;
    sn =  0; sr1 =  0; sr2 =  0;
    deltan = 0; epsilonn = 0; 
    gamma = 0; gammal = 0; gammal2 = 0;
    eta = 0; etal = 0; etal2 = 0;
    vepln = 0; veplnl = 0; veplnl2 = 0; ul3= 0;
    ul2 = 0; ul = 0; u = 0; rnorm = betan;
    xl2norm = 0; 
    Anorm = 0; Acond = 1;    
    relres = 1;
    QLPiter = 0;
    
    x.setZero(n,m); w.setZero(n,m); wl.setZero(n,m); 
    resvec(0,0) = beta1;

    //NB User input
    maxxnorm = 1e7;

    
    m_error = r2.norm() / m_rhs_norm;
    if (m_error < m_tol)
        return true;

    return false;
}

template<class T>
bool gsMinResQLP<T>::step( typename gsMinResQLP<T>::VectorType& x )
{
    betal = beta;
    beta= betan;
    v = (1.0/beta)*r3;
    m_mat->apply(v,r3);
    //Assume that m_num_iter starts at 1 (at this point in the code)
    //TODO check assumption 
    if (m_num_iter > 1)
        r3 = r3 - (beta/betal)*r1;
    
    alpha = T (r3.col(0).dot(v.col(0)));
    r3 = r3 - (alpha/beta)*r2;
    r1.swap(r2);
    r2.swap(r3);
    
    m_precond->apply(r2,r3);
    betan = r2.col(0).dot(r3.col(0));
    GISMO_ASSERT(betan > 0, "Preconditioner is indefinite");
    betan = math::sqrt(betan);
    
    //pnorm = rho
    pnorm = math::sqrt(alpha*alpha + betal*betal + betan*betan); 
    
    dbar = deltan;
    delta = cs*dbar + sn*alpha;
    epsilon = epsilonn;
    gbar = sn*dbar - cs*alpha; 
    epsilonn= sn*betan;
    deltan = -cs*betan; 

    gammal3 = gammal2;
    gammal2 = gammal;
    gammal = gamma; 
    MinResQLPSymOrtho(gbar, betan, cs, sn, gamma);
    //real_t gammaTmp = gamma; used for (normal) minres update
    taul2 = taul;
    taul  = tau;
    tau   = cs*phi; 
    Axnorm =  math::sqrt(Axnorm*Axnorm + tau*tau); 
    phi   = sn*phi;

    
    if (m_num_iter > 2)
    {
        veplnl2 = veplnl;
        etal2 = etal; 
        etal = eta;
        real_t delta_tmp = sr2*vepln - cr2*delta;
        veplnl = cr2*vepln + sr2*delta;
        delta = delta_tmp;
        eta = sr2*gamma;
        gamma = -cr2*gamma;
    }
    if (m_num_iter > 1)
    {
        MinResQLPSymOrtho(gammal, delta, cr1, sr1, gammal);
        vepln = sr1*gamma;
        gamma = -cr1*gamma;
    }

    
    real_t ul4 = ul3;
    ul3 = ul2;
    if (m_num_iter > 2)
    {
        ul2 = (taul2 - etal2*ul4 - veplnl2*ul3) / gammal2;
    }
    if (m_num_iter > 1)
    {
        ul = ( taul  - etal *ul3 - veplnl *ul2) / gammal;    
    }
    real_t xnorm_tmp = math::sqrt(xl2norm*xl2norm + ul2*ul2 + ul*ul); 
    if ( math::abs(gamma) > std::numeric_limits<real_t>::min() && xnorm_tmp < maxxnorm)
    {
        u = (tau - eta*ul2 - vepln*ul) / gamma; 
        if ( math::sqrt(xnorm_tmp*xnorm_tmp + u*u) > maxxnorm)
        {
            u = 0; //NB flag = 6 LINE 374
        }        
    }
    else
    {
        u = 0; //NB flag = 9 LINE 379
    }
    xl2norm = math::sqrt(xl2norm*xl2norm + ul2*ul2);
    
    xnorm   = math::sqrt(xl2norm*xl2norm + ul*ul + u*u);
    
    //MINRES-QLP updates 
    //TODO: add IF/ELSE if we want to integrate
    //      normal MINRES with MINRES-QLP
    QLPiter += 1;
    if (QLPiter == 1)
    {
        xl2.setZero(m_mat->cols(),1);
        if (m_num_iter > 1)
        {
            if (m_num_iter > 3)
            {
                wl2 = gammal3*wl2 + veplnl2*wl + etal*w;
            }
            if (m_num_iter > 2)
            {
                wl = gammal_QLP*wl + vepln_QLP*w;
            }
            w = gamma_QLP*w;
            xl2 = x - wl*ul_QLP - w*u_QLP;
        }
    }
    if (m_num_iter == 1)
    {
        wl2.swap(wl);
        wl = v*sr1;
        w  = -v*cr1;
    }
    else if (m_num_iter ==2)
    {
        wl2.swap(wl);
        wl  = w*cr1 + v*sr1;
        w   = w*sr1 - v*cr1;
    }
    else
    {
        wl2.swap(wl);   
        wl.swap(w);
        w  = wl2*sr2 - v*cr2;
        wl2 = wl2*cr2 + v*sr2;
        v  = wl *cr1 + w*sr1;    
        w   = wl *sr1 - w*cr1;
        wl.swap(v);
    }
    xl2 = xl2 + wl2*ul2;
    x   = xl2 + wl *ul + w*u;

    real_t gammal_tmp = gammal;
    MinResQLPSymOrtho(gammal, epsilonn, cr2, sr2, gammal);

    gammal_QLP = gammal_tmp;     
    vepln_QLP = vepln;     
    gamma_QLP = gamma; 
    ul_QLP    = ul;
    u_QLP     = u;

    abs_gamma = math::abs(gamma);
    
    real_t tmp_max  = std::max(gammal, abs_gamma);
    Anorm = std::max(Anorm, pnorm);
    Anorm = std::max(Anorm, tmp_max);
    if (m_num_iter == 1)
    {
        gmin = gamma; 
        gminl = gmin;
    }
    else if (m_num_iter > 1)
    {
        gminl2 = gminl;
        gminl = gmin;
        
        gmin = std::min(gminl2, gammal);
        gmin = std::min(gmin, abs_gamma);
    }

    
    Acond   = Anorm/gmin;
    rnorml  = rnorm;
    if (phi != 0.0) //LINE 441
        rnorm = phi; //NB if flag != 9;
    relres   = rnorm / (Anorm*xnorm + beta1);
    rootl    = math::sqrt(gbar*gbar + deltan*deltan);   
    Arnorml  = rnorml*rootl;
    relAresl = rootl / Anorm; 

    
    // Test for convergence
    if (m_inexact_residual)
    {
        m_error = relres;
    }
    else
    {
        m_mat->apply(x,negResidual);
        gsMatrix<T> r1_tmp = m_rhs - negResidual;
        m_error = r1_tmp.norm()/m_rhs_norm;
    }
    if (m_error < m_tol)
        return true;

    
    resvec(m_num_iter, 0) = rnorm;
    Aresvec(m_num_iter-1, 0)  = Arnorml;

    return false;

}

template<class T>
void gsMinResQLP<T>::finalizeIteration( typename gsMinResQLP<T>::VectorType& x)
{
    m_mat->apply(x,negResidual);
    r1 = m_rhs - negResidual;
    rnorm = r1.norm();
    m_mat->apply(r1, negResidual);
    Arnorm = negResidual.norm();
    xnorm = x.norm();
    relres = rnorm/(Anorm*xnorm + beta1);
    m_error = relres;
    relAres = 0; 
    if(rnorm > std::numeric_limits<real_t>::min())
    {
        relAres = Arnorm / (Anorm*rnorm);
    }
    Aresvec(m_num_iter, 0)  = Arnorm;
    Aresvec.conservativeResize(m_num_iter+1,1);
    resvec.conservativeResize(m_num_iter+1,1);
    r1.clear(); r2.clear(); r3.clear(); 
    
    negResidual.clear(); v.clear();
    xl2.clear();
    wl.clear(); w.clear(); 
    wl2.clear(); 
}


}
