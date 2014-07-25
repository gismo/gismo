
#pragma once 

#include <gsUtils/gsDeCasteljau.hpp>
#include <gsUtils/gsSparseRows.hpp>
#include <gsDataStructures/gsMesh/gsMesh.h>

namespace gismo
{

template<class T>
gsBernsteinBasis<T>& gsBernsteinBasis<T>::component(unsigned i) const
{
    if ( i == 0 )
        return const_cast<gsBernsteinBasis&>(*this);
    else
        GISMO_ERROR("gsBernsteinBasis has only one component");
}

template <class T>
void gsBernsteinBasis<T>::connectivity(const gsMatrix<T> & nodes, 
                                       gsMesh<T> & mesh) const
{
    const index_t sz  = size();
    GISMO_ASSERT( nodes.rows() == sz, "Invalid input.");

    // Add adges and vertices
    for(index_t i = 0; i< sz - 1; ++i )
        mesh.addEdge( nodes.row(i  ).transpose(),
                      nodes.row(i+1).transpose());
}


template <class T>
void gsBernsteinBasis<T>::anchors_into(gsMatrix<T> & result) const
{
    result.resize(1, this->size() );
    int s = 0;
    for ( int i = 1; i< m_breaks.size(); ++i)
        for ( int j = 0; j< m_p; ++j)
            result(0, s++) = m_breaks[i-1] + 
                ( m_breaks[i] - m_breaks[i-1] ) * T(j)/T(m_p) ;

    result(0, s) = m_breaks.last();
}


template<class T>
void gsBernsteinBasis<T>::active_into(const gsMatrix<T>& u, gsMatrix<unsigned>& result) const 
{
    result.resize(m_p+1, u.cols());
    for (index_t j = 0; j < u.cols(); ++j)
    {
        unsigned first = firstActive(u(0,j));
        for (int i = 0; i <= m_p; ++i)
            result(i,j) = first + i;
    }
}


template<class T> 
gsMatrix<unsigned> * gsBernsteinBasis<T>::boundary() const
{
    gsMatrix<unsigned> * res = new gsMatrix<unsigned>(2,1);
    (*res)(0,0)= 0;
    (*res)(1,0)= this->size() -1;
    return res;
}


template<class T> 
gsMatrix<unsigned> * gsBernsteinBasis<T>::boundary(boundary::side const & s ) const
{
    gsMatrix<unsigned> * res = new gsMatrix<unsigned>(1,1);
    switch (s) {
    case boundary::left : // left
        (*res)(0,0)= 0;
        break;
    case boundary::right : // right
        (*res)(0,0)= this->size() -1;
        break;       
    default:
        GISMO_ERROR("gsBernsteinBasis: valid sides is left or right.");
    };            
    return res;
}


template<class T>
gsMatrix<T> gsBernsteinBasis<T>::support() const 
{
    gsMatrix<T> res(1,2);
    res << domainStart() , domainEnd() ;
    return res ;
}

template<class T>
gsMatrix<T> gsBernsteinBasis<T>::support(const unsigned & i) const 
{
    gsMatrix<T> res(1,2);
    res << m_breaks[ i % m_p ] , m_breaks[ i % m_p + 1 ] ;
    return res ;
}


template<class T> 
void gsBernsteinBasis<T>::eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const 
{
    result.resize(m_p+1, u.cols());

    for ( index_t k=0; k< u.cols(); ++k )
    {
        const int i = m_breaks.findspan( u(0,k) ) ;
        // Check if point is valid
        // assert( u(0,k) >= m_breaks.first() && u(0,k) <= m_breaks.last() ) ;
      
        result(0,k)= T(1);
        const T u0 = ( u(0,k) - m_breaks[i])/(m_breaks[i+1]-m_breaks[i]) ;
        const T u1 = T(1.0) - u0 ;
        for (int j=1;j<=m_p; ++j)
        {
            T saved = 0.0;
            for (int r=0;r<j; ++r)
            {
                const T temp = result(r,k);
                result(r,k)= saved + u1*temp;
                saved= u0 * temp ;
            }
            result(j,k)= saved;
        }
    }
}
    

template<class T> 
void gsBernsteinBasis<T>::evalSingle_into(unsigned i, const gsMatrix<T> & u, gsMatrix<T>& result) const {

    result.resize(1, u.cols() );
    STACK_ARRAY(T, N, m_p + 1);

    for (index_t s=0;s<u.cols(); ++s)
    {
        // Locality property
        unsigned fs = m_p * m_breaks.findspan(u(0,s)) ;
        if (   i < fs || i > fs +m_p )
        {
            result(0,s)= T(0);
            continue;
        }

        // Initialize zeroth degree functions
        for (int j=0;j<=m_p; ++j)
            N[j] = 0;
        N[m_p-i] = T(1);
        T u0 = ( u(0,s) - m_breaks[i])/(m_breaks[i+1]-m_breaks[i]) ;
        T u1 = T(1.0) - u0 ;
        for (int k=1;k<=m_p; ++k)
            for (int j=m_p;j>=k; --j)
                N[j]= u1  * N[j] + N[j-1] * u0 ;
        result(0,s)= N[m_p];
    }
}


template<class T> inline
void gsBernsteinBasis<T>::eval_into(const gsMatrix<T> &u, const gsMatrix<T> & coefs, gsMatrix<T>& result) const 
{  
    gsDeCasteljauEval(u,m_breaks, m_p, coefs, result);

    //gsHornerBezier(u, m_breaks, m_p, coefs, result);// Has a bug ( check partition of unity)

    //gsBasis<T>::eval_into(u,coefs,result); // defalult gsBasis implementation
}


template<class T> inline
void gsBernsteinBasis<T>::deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const 
{
    assert( u.rows() == 1 ) ;
    gsMatrix<T> tmp;
    this->evalAllDers_into(u,1,tmp);
    result = tmp.block(m_p+1,0,m_p+1,u.cols());
}

template<class T>  inline
void gsBernsteinBasis<T>::derivSingle_into(unsigned & i, const gsMatrix<T> & u, gsMatrix<T>& result ) const 
{
    // TO DO
    gsBasis<T>::derivSingle_into(i, u, result);
}


template<class T> inline 
void gsBernsteinBasis<T>::deriv_into(const gsMatrix<T> & u, const gsMatrix<T> & coefs, gsMatrix<T>& result ) const 
{ 
    // TO DO specialized computation for gsBernsteinBasis
    gsBasis<T>::deriv_into(u,coefs,result);
}


template<class T> inline
gsMatrix<T> * gsBernsteinBasis<T>::deriv2(const gsMatrix<T> & u ) const 
{
    return new gsMatrix<T>(this->evalAllDers(u,2)->block(2*m_p+2,0,m_p+1,u.cols()) );
}


template<class T> inline
gsMatrix<T> * gsBernsteinBasis<T>::laplacian(const gsMatrix<T> & u ) const 
{
    gsMatrix<T> * der2 = new gsMatrix<T>(this->evalAllDers(u,2)->block(2*m_p+2,0,m_p+1,u.cols()) ); 
    gsMatrix<T> * ev = new gsMatrix<T>( der2->colwise().sum() ) ;
    delete der2;
    return ev;
}

template<class T>
gsBernsteinBasis<T> * gsBernsteinBasis<T>::clone() const
{ return new gsBernsteinBasis(*this); };


template<class T>
void gsBernsteinBasis<T>::evalAllDers_into(const gsMatrix<T> & u, int n, gsMatrix<T>& result) const  
{
    assert( u.rows() == 1 );
    const int p1 = m_p + 1;       // degree plus one

    STACK_ARRAY(T, ndu,  p1 * p1 );
    STACK_ARRAY(T, a, 2 * p1);

    result.resize( (n+1)*(m_p + 1), u.cols() ) ;  

    for (index_t v = 0; v < u.cols(); ++v) // for all columns of u
    {
        // Run evaluation algorithm but keep the triangle
        unsigned i = m_breaks.findspan( u(0,v) ) ;
    
        ndu[0] = T(1) ;
        T u0 = ( u(0,v) - m_breaks[i])/ (m_breaks[i+1]-m_breaks[i]) ;
        T u1 = T(1.0) - u0 ;
        for(int j=1; j<= m_p; j++)
        {
            T saved = T(0) ;
            for(int r=0; r<j ; r++)
            {
                // Lower triangle
                ndu[j*p1 + r] = u1 + u0 ;
                const T temp = ndu[r*p1 + j-1] / ndu[j*p1 + r] ;
                // Upper triangle
                ndu[r*p1 + j] = saved + u1 * temp ;
                saved = u0 * temp ;
            }  
            ndu[j*p1 + j] = saved ;
        }
      
        // Assign 0-derivative equal to function values
        //result.block(0,v, p1,1) = ndu.col(m_p);
        for (int j=0; j <= m_p ; ++j )
            result(j,v) = ndu[j*p1 + m_p];
    
        // Compute the derivatives
        for(int r = 0; r <= m_p; r++)
        {
            // alternate rows in array a
            T* a1 = &a[0];
            T* a2 = &a[p1];

            a1[0] = T(1) ;

            // Compute the k-th derivative of the r-th basis function
            T tmp = (m_breaks[i+1]-m_breaks[i]);
            for(int k=1; k<=n; k++)
            {
                int rk,pk,j1,j2 ;
                T d(0) ;
                rk = r-k ; pk = m_p-k ;
        
                if(r >= k)
                {
                    a2[0] = a1[0] / ndu[ (pk+1)*p1 + rk] ;
                    d = a2[0] * ndu[rk*p1 + pk] ;
                }
        
                if (rk >= -1)  j1 = 1;   else j1 = -rk;
                if (r-1 <= pk) j2 = k-1; else j2 = m_p - r;
        
                for(int j = j1; j <= j2; j++)
                {
                    a2[j] = (a1[j] - a1[j-1]) / ndu[(pk+1)*p1 + rk+j] ;
                    d += a2[j] * ndu[(rk+j)*p1 + pk] ;
                }
        
                if(r <= pk)
                {
                    a2[k] = -a1[k-1] / ndu[(pk+1)*p1 + r] ;
                    d += a2[k] * ndu[r*p1 + pk] ;
                }
                result(k*p1 + r, v) = d / tmp ;
                tmp *= tmp;
                std::swap(a1, a2);              // Switch rows
            }
        }
    }// end v

    // Multiply through by the correct factors
    int r = m_p ;
    for(int k=1;k<=n;k++)
    {
        result.middleRows(k*p1, p1).array() *= T(r) ;
        r *= m_p - k ;
    }
}

template <class T>
void gsBernsteinBasis<T>::uniformRefine_withCoefs(gsMatrix<T> & coefs, int numKnots)
{
    std::vector<T> u = this->m_breaks.get();

    for (unsigned i = 0; i < u.size() - 1; ++i)
        for (int k = 1; k <= numKnots; ++k)
            gsDeCasteljauSubdivide(m_breaks, m_p, coefs, ((numKnots+1-k) * u[i] + k * u[i+1]) / (numKnots + 1));
}

template <class T>
void gsBernsteinBasis<T>::uniformRefine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, int numKnots)
{
    std::vector<T> u = this->m_breaks.get() ; 
    gsSparseRows<T> trans;
    trans.setIdentity( this->size() );
    
    for (unsigned i = 0; i < u.size() - 1; ++i)
        for (int k = 1; k <= numKnots; ++k)
            gsDeCasteljauSubdivide(m_breaks, m_p, trans, ((numKnots+1-k) * u[i] + k * u[i+1]) / (numKnots + 1));

    trans.toSparseMatrix( transfer );
}


template <class T>
void gsBernsteinBasis<T>::insertKnot_withCoefs(T const knot, gsMatrix<T> & coefs)
{ 
    gsDeCasteljauSubdivide(m_breaks, m_p, coefs, knot); 
}



} // namespace gismo
