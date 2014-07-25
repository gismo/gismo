
#pragma once

#include <gsUtils/gsDebug.h>
# include <gsNurbs/gsBSpline.h>

#include <numeric>

namespace gismo {


template<class T>
template<class KnotVectorType>
void gsBSplineSolver<T>::initSolver(gsBSpline<T,KnotVectorType> const & bsp , 
                                    int const & coord, T const & tr, 
                                    T const & tol, unsigned const & N)
{
    m_n  = bsp.coefsSize();
    m_d  = bsp.degree();
    m_k  = 1;
    eps  = tol;
    x    = 0.0;
    maxn = N + m_n;
  
    const KnotVectorType & kv = bsp.knots();
  
    m_t.resize(m_n+N+m_d+1); 
    m_c.resize(m_n+N);  

    for ( unsigned i = 0; i < m_n; ++i ) 
    {
        m_c[i]= bsp.coef(i, coord) - tr;
        m_t[i]= kv[i];
    }
  
    for ( unsigned i = m_n; i < m_n + m_d +1 ; ++i ) 
        m_t[i]= kv[i];
}
    
template<class T>
bool gsBSplineSolver<T>::firstRoot (gsBSpline<T,gsKnotVector<T> > const & bsp , 
                                    int const & coord, T const & tr, 
                                    T const & tol, unsigned const & N)
{
    initSolver<gsKnotVector<T> >(bsp,coord,tr,tol,N);
    return nextRoot();
}

template<class T>
bool gsBSplineSolver<T>::firstRoot (gsBSpline<T, gsCompactKnotVector<T> > const & bsp,
                                    int const & coord, T const & tr, T const & tol, 
                                    unsigned const & N)
{
    initSolver<gsCompactKnotVector<T> >(bsp,coord,tr,tol,N);
    return nextRoot();
}

template<class T>
void gsBSplineSolver<T>::allRoots (gsBSpline<T,gsKnotVector<T> > const & bsp, 
                                   std::vector<T> & result, 
                                   int const & coord, T const & tr, T const & tol, 
                                   unsigned const & N)
{
    result.clear();
    initSolver<gsKnotVector<T> > (bsp,coord,tr,tol,N);

    while ( nextRoot() )
    {
        result.push_back( value() ) ;
    }
    //gsLog<< "gsBSplineSolver: found "<< result.size() <<" roots.\n  ";
}

template<class T>
void gsBSplineSolver<T>::allRoots (gsBSpline<T,gsCompactKnotVector<T> > const & bsp, 
                                   std::vector<T> & result, 
                                   int const & coord, T const & tr, T const & tol, 
                                   unsigned const & N)
{
    result.clear();
    initSolver<gsCompactKnotVector<T> > (bsp,coord,tr,tol,N);

    while ( nextRoot() )
    {
        result.push_back( value() ) ;
    }
    //gsLog<< "gsBSplineSolver: found "<< result.size() <<" roots.\n  ";
}

template<class T>
int gsBSplineSolver<T>::insertKnot(int mu) 
{
    // find span (in case mu is not tight)
    mu = math::max(mu,m_d);

    while ( x>=m_t[mu+1]) mu++; 
  
    // make room for one coefficient at position mu
    m_c.insert(m_c.begin()+mu,0);
    //for ( int i=m_n; i>mu; i--) m_c[i] = m_c[i-1];

    // Compute new coefficients
    for ( int i=mu; i>=mu-m_d+1; i--) 
    {
        const T alpha = (x-m_t[i])/(m_t[i+m_d]-m_t[i]);
        m_c[i] = (1.0-alpha) * m_c[i-1] + alpha * m_c[i];
    }

    // set last knot
    m_t[m_n+m_d] = m_t[m_n+m_d-1];

    //insert the knot in the knot sequence
    m_t.insert(m_t.begin()+mu+1,x);
    // for ( int i=m_n; i>mu+1; --i)
    //     m_t[i] = m_t[i-1];
    // m_t[mu+1] = x;
  
    m_n++;
    return mu;
}
  
template<class T>
bool gsBSplineSolver<T>::nextRoot()
{
    // Do until maxn knots are inserted
    while  ( m_n < maxn ) 
    {
        // Find sign change
        while( m_k<m_n && (m_c[m_k-1]*m_c[m_k] > 0) )
        {
            ++m_k;
            //gsLog<<"m_c[m_k]= "<< m_c[m_k] <<" ( m_k="<<m_k<<")\n";    
        }

        // No sign change ?
        if ( m_k >= m_n )
            return false;
        //gsLog<<"Sign change at "<< m_c[m_k]<<", "<< m_c[m_k+1] <<" (position "<<m_k<<")\n";

        // Interval converged ?
        const T diff = m_t[m_k+m_d]-m_t[m_k];
        if (diff<eps) 
        {
            x = m_t[m_k];// root found
            //gsLog<<"diff converged: Root found at "<< m_k<<", val="<< x <<"m_n="<<m_n <<"\n";
            break;
        }

        // Horizontal alignment ?
        const T cdiff = m_c[m_k-1]-m_c[m_k];
        if (math::abs(cdiff)<eps) 
        {
            ++m_k;
            continue;
        }

        // compute line slope
        const T lambda = m_c[m_k-1] / cdiff ;

        // compute intersection 
        x = (std::accumulate(m_t.begin()+m_k, m_t.begin()+m_k+m_d, T(0.0) )
             + lambda*diff ) / T(m_d);
        //x = ( m_t.segment(m_k, m_d).sum() + lambda*diff ) / T(m_d);

        // Stopping criterion
        const T e = math::max(x, m_t[m_k+m_d-1]) - math::min(x, m_t[m_k+1] );
        if (e < eps)
        {
            //gsLog<<"eps: Root found after knot "<< m_k <<" is "<<m_t[m_k]<<", val="<< x <<"\n";
            break;
        }

        insertKnot(m_k);
        //gsLog<<"Inserted knot, m_k="<< m_k <<", m_t[m_k]="<<m_t[m_k]<<", m_n="<< m_n<<"\n";  
    
    } // end for
  

    //gsLog<<"Loop finished with "<< x <<" m_k="<<m_k<<", m_c[m_k]="<<m_c[m_k]<<", m_t[m_k]="<<m_t[m_k]<< ", m_n="<< m_n<<"\n";

    if ( m_n>=maxn )
    {
        gsWarn<<"gsBSplineRoot: Maximum number of knots reached: "<< m_n << "\n";
        //gsInfo<< "m_c:\n"<< m_c.transpose()<<"\n";
        //gsInfo<< "m_t:\n"<< m_t.transpose()<<"\n";
        return false;
    }

    m_k++; // next root

    while ( m_k<m_n && math::fabs( m_t[m_k]-x)<eps ) m_k++; 

    //gsLog<<"next m_k="<< m_k <<" is "<< m_t[m_k]<<"\n";
    return true;
}


}; //namespace gismo
