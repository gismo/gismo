/** @file gsTensorBasis.hpp

    @brief Provides implementation of TensorBasis class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsTensor/gsTensorTools.h>

#include <gsTensor/gsTensorDomainIterator.h>
#include <gsTensor/gsTensorDomainBoundaryIterator.h>

#include <gsCore/gsBoundary.h>
#include <gsUtils/gsMesh/gsMesh.h>

#if defined __INTEL_COMPILER
#pragma warning (disable : 175 ) /* disable subscript out of range warning */
#endif


namespace gismo
{

template<unsigned d, class Basis_t >
gsTensorBasis<d,Basis_t>::gsTensorBasis( Basis_t* x,  Basis_t* y) 
{ 
    GISMO_ASSERT( d==2, "gsTensorBasis error: wrong dimension." );

    if ( x->dim() == unsigned(1) && y->dim() == unsigned(1) )
    {
        m_bases[0] = x;
        m_bases[1] = y;
    } else 
        GISMO_ERROR("gsTensorBasis error: Spaces must be of topological dimension 1.");
}


template<unsigned d, class Basis_t >
gsTensorBasis<d,Basis_t>::gsTensorBasis(Basis_t* x, Basis_t* y, Basis_t* z) 
{ 
    GISMO_ASSERT( d==3, "gsTensorBasis error: wrong dimension." );
  
    GISMO_ASSERT( x->dim() == 1 && y->dim() == 1 && z->dim() == 1,
                  "Spaces must be of topological dimension 1." );
  
    if ( d==3 )
    {
        m_bases[0] = x;
        m_bases[1] = y;
        m_bases[2] = z;
    } 
    else 
        GISMO_ERROR("gsTensorBasis incorrect constructor for "<<d<<"D basis");
}

template<unsigned d, class Basis_t >
gsTensorBasis<d,Basis_t>::gsTensorBasis(Basis_t* x, Basis_t* y, Basis_t* z, Basis_t* w) 
{ 
    GISMO_ASSERT( d==4, "gsTensorBasis error: wrong dimension." );
  
    GISMO_ASSERT( x->dim() == 1 && y->dim() == 1 && z->dim() == 1 && w->dim() == 1,
                  "Spaces must be of topological dimension 1." );
  
    if ( d==4 )
    {
        m_bases[0] = x;
        m_bases[1] = y;
        m_bases[2] = z;
        m_bases[3] = w;
    } 
    else 
        GISMO_ERROR("gsTensorBasis incorrect constructor for "<<d<<"D basis");
}


template<unsigned d, class Basis_t >
gsTensorBasis<d,Basis_t>::gsTensorBasis( std::vector<Basis_t*> const & bb )
{ 
    GISMO_ASSERT( d==bb.size(), "gsTensorBasis error: wrong number of bases given ("
                  << bb.size()<< ", expected "<< d );

    for (unsigned i = 0; i < d; ++i)
        m_bases[i] = bb[i];
}


/// Copy Constructor
template<unsigned d, class Basis_t >
gsTensorBasis<d,Basis_t>::gsTensorBasis( const gsTensorBasis & o)
: gsBasis<T>(o)
{
    for (unsigned i = 0; i < d; ++i)
        m_bases[i] = o.m_bases[i]->clone();
}


template<unsigned d, class Basis_t >
gsTensorBasis<d,Basis_t>& gsTensorBasis<d,Basis_t>::operator=( const gsTensorBasis & o)
{
    if ( this == &o )
        return *this;

    gsBasis<T>::operator=(o);

    for (unsigned i = 0; i < d; ++i)
    {
        delete m_bases[i];
        m_bases[i] = o.m_bases[i]->clone();
    }
    return *this;
}


template<unsigned d, class Basis_t >
void gsTensorBasis<d,Basis_t>::anchors_into(gsMatrix<T>& result) const
{
    gsMatrix<T> gr[d];
    gsVector<unsigned, d> v, size;
    result.resize( d, this->size() );

    for (unsigned i = 0; i < d; ++i)
    {
        gr[i] = m_bases[i]->anchors();
        size[i] = this->size(i);
    }

    // iterate over all tensor product basis functions
    v.setZero();
    unsigned r = 0;
    do {
        // Make tensor product of greville points
        for (unsigned i=0; i<d; ++i )
            result(i,r)=  (gr[i])( 0, v(i) );
        ++r ;
    } while (nextLexicographic(v, size));
}

template<unsigned d, class Basis_t >
void gsTensorBasis<d,Basis_t>::connectivity(const gsMatrix<T> & nodes, 
                                            gsMesh<T> & mesh) const
{
    const index_t sz  = size();
    GISMO_ASSERT( nodes.rows() == sz, "Invalid input.");

    // Add vertices
    for(index_t i = 0; i< sz; ++i )
        mesh.addVertex( nodes.row(i).transpose() );

    // Starting from lower corner
    const gsVector<unsigned,d> low = gsVector<unsigned,d>::Zero();    

    // Last tensor-index
    gsVector<unsigned, d> end;
    for (unsigned i = 0; i < d; ++i)
        end(i) = this->size(i) - 1;

    unsigned k, s;
    gsVector<unsigned,d> v, upp;
    for (unsigned i = 0; i < d; ++i) // For all axes
    {
        s      = stride(i);
        v      = low;
        upp    = end;
        upp[i] = 0; // suppress to face v[i]==0

        do // Insert all edges normal to axis i 
        {
            k = index(v);
            for (unsigned j = 0; j != end[i]; ++j)
            {
                mesh.addEdge(k, k + s );
                k += s;
            }
        }
        while ( nextCubePoint(v, low, upp) );
    }
}

// generic version for the case that the basis does not implement firstActive() / numActive():
//
//~ template<unsigned d, class Basis_t >
//~ gsMatrix<unsigned> * gsTensorBasis<d,Basis_t>::active(const gsMatrix<T> & u ) const
//~ {
//~   gsMatrix<unsigned>* act[d];
//~   gsVector<unsigned, d> v, size;
//~ 
//~   // Get univariate active basis functions
//~   unsigned nb = 1;
//~   for (unsigned i = 0; i < d; ++i)
//~   {
//~     gsMatrix<unsigned>* curAct = m_bases[i]->active(u.row(i)) ;
//~     act[i] = curAct;
//~     nb *= curAct->rows();
//~     size[i]= curAct->rows();
//~   }
//~ 
//~   gsMatrix<unsigned> * res= new gsMatrix<unsigned>( nb, u.cols() );  
//~   
//~   // iterate over all tensor product active functions
//~   unsigned r = 0;
//~   v.setZero();
//~   do {
//~     // Fill with active bases indices
//~     for ( index_t j=0; j<u.cols(); ++j)
//~     {
//~       nb= (*act[d-1])( v(d-1) ,j) ;//compute global index in the tensor product
//~       for ( int i=d-2; i>=0; --i )
//~         nb = nb * m_bases[i]->size() + (*act[i])( v(i),j) ;
//~       (*res)( r, j )= nb;
//~     }
//~     ++r ;
//~   } while (nextLexicographic(v, size));
//~   
//~   // free temporary memory
//~   for (unsigned i = 0; i < d; ++i)
//~     delete act[i];
//~ 
//~   return res;
//~ };

/*
 * NB:
 * This implementation assumes that the component bases implement the firstActive() / numActive() protocol.
 * In particular, their active basis functions must always be continuous intervals.
 * This is the case for all current component bases, so we only keep this version for now.
 * Above, commented out, is the generic version which is quite a bit slower.
 */
template<unsigned d, class Basis_t >
void gsTensorBasis<d,Basis_t>::active_into(const gsMatrix<T> & u, gsMatrix<unsigned>& result) const
{
    unsigned firstAct[d];
    gsVector<unsigned, d> v, size;

    // count active functions in each tensor direction
    unsigned numAct = 1;
    for (unsigned i = 0; i < d; ++i)
    {
        size[i] = m_bases[i]->numActive();
        numAct *= size[i];
    }

    result.resize( numAct, u.cols() );  
  
    // Fill with active bases indices
    for (index_t j = 0; j < u.cols(); ++j)
    {
        // get the active basis indices for the component bases at u(:,j)
        for (unsigned i = 0; i < d; ++i)
        {
            firstAct[i] = m_bases[i]->firstActive( u(i,j) );
        }

        // iterate over all tensor product active functions
        unsigned r = 0;
        v.setZero();
        do
        {
            int gidx = firstAct[d-1] + v(d-1);    //compute global index in the tensor product
            for ( int i=d-2; i>=0; --i )
                gidx = gidx * this->size(i) + firstAct[i] + v(i);

            result(r, j) = gidx;
            ++r ;
        } while (nextLexicographic(v, size));
    }
}

template<unsigned d, class Basis_t >
void gsTensorBasis<d,Basis_t>::active_cwise(
    const gsMatrix<typename gsTensorBasis<d,Basis_t>::T> & u, 
    gsVector<unsigned,d>& low, 
    gsVector<unsigned,d>& upp ) const
{
    for (index_t j = 0; j < u.cols(); ++j)
    {
        for (unsigned i = 0; i < d; ++i)
        {
            low[i] = m_bases[i]->firstActive( u(i,j) );
            upp[i] = low[i] + m_bases[i]->degree();
        }
    }
}


template<unsigned d, class Basis_t >
typename gsMatrix<unsigned>::uPtr gsTensorBasis<d,Basis_t>::slice(int dir, int k) const
{
    GISMO_ASSERT( dir>=0 &&  dir < this->dim(), "Invalid slice direction requested" );
    GISMO_ASSERT( k >=0 &&  k < this->trueSize(dir), "Invalid slice position requested" );

    unsigned sliceSize = 1, r = 0;
    gsVector<unsigned, d> low, upp;

    // Set low and upp to point to the range of indices
    for (unsigned i = 0; i < d; ++i)
    {
        sliceSize *= this->size(i); // or trueSize(i) ?
        low(r) = 0;
        upp(r++) = this->size(i); // or trueSize(i) ?
    }
    sliceSize /= upp(dir);
    upp(dir) = k + 1;
    low(dir) = k;

    typename gsMatrix<unsigned>::uPtr res ( new gsMatrix<unsigned>(sliceSize,1) );

    // iterate over all tensor product basis indices
    gsVector<unsigned,d> v = low;
    r = 0;
    do {
        (*res)(r++,0) = this->index(v);
    } while ( nextLexicographic(v, low, upp) );

    return res;
}


template<unsigned d, class Basis_t >
gsMatrix<unsigned> * gsTensorBasis<d,Basis_t>::boundary() const
{
    gsMatrix<unsigned> bd;
    std::set<unsigned> bdofs;

    for (unsigned k = 0; k != d; ++k)
    {
        bd = this->slice(k, 0);
        for (index_t i = 0; i < bd.size(); ++i)
            bdofs.insert( bd(i) );

        bd = this->slice(k, size(k) - 1);
        for (index_t i = 0; i < bd.size(); ++i)
            bdofs.insert( bd(i) );
    }

    return makeMatrix<unsigned>(bdofs.begin(), bdofs.size(), 1 ).release();

    /* // returns boundary with repetitions
       unsigned sz(0), i(0), r(0);
       unsigned bdofs, dofs = this->size() ;

       //compute size of boundary (with repetitions at the corners)
       for (unsigned k = 0; k < d; ++k)
       sz += 2 * dofs / size(k);
       gsMatrix<unsigned> * res = new gsMatrix<unsigned>(sz,1);

       // Fill boundary DoFs by taking slices
       r= sz ;
       for (unsigned k = 0; k < d; ++k)
       {
       bdofs= dofs / size(k);
       r-=  bdofs;
       res->block(r,0, bdofs, 1) = *safe( this->slice(i, size(k) - 1) );
       r-= bdofs ;
       res->block(r,0, bdofs, 1) = *safe( this->slice(i,0) );
       i++;
       }

       return res;
    */
};


template<unsigned d, class Basis_t >
gsMatrix<unsigned> * gsTensorBasis<d,Basis_t>::boundary(boxSide const& s) const
{
    //get m_bases index and start or end case
    int k = direction(s);
    int r = parameter(s);
    return this->slice(k, (r ? size(k) - 1 : 0) ).release();
};

template <unsigned d, class BB, class B>
struct MakeBoundaryBasis
{
    static BB* make (const std::vector< B * >& bases)
    {
        return new BB( bases );
    }
};

template <class BB, class B>
struct MakeBoundaryBasis<2, BB, B>
{
    static BB* make (const std::vector< B * >& bases)
    {
        return bases[0];
    }
};


template<unsigned d, class Basis_t >
typename gsTensorBasis<d,Basis_t>::BoundaryBasisType * 
gsTensorBasis<d,Basis_t>::boundaryBasis(boxSide const& s) const
{   
    unsigned dir = direction( s );

    std::vector<Basis_t*> rr;
    rr.reserve( d - 1 );
    for ( unsigned i=0; i < d; ++i )
        if (i != dir)
            rr.push_back( m_bases[i]->clone() );
    
    return MakeBoundaryBasis<d, BoundaryBasisType, Basis_t>::make( rr );
}


template<unsigned d, class Basis_t >
gsMatrix<typename gsTensorBasis<d,Basis_t>::T> gsTensorBasis<d,Basis_t>::support() const 
{
    gsMatrix<T> res(d,2);
    for (unsigned i = 0; i < d; ++i)
        res.row(i) =  m_bases[i]->support();
    return res;
}

template<unsigned d, class Basis_t >
gsMatrix<typename gsTensorBasis<d,Basis_t>::T> 
gsTensorBasis<d,Basis_t>::support(const unsigned & i) const 
{
    gsMatrix<T> res(d,2);
    gsVector<unsigned, d> ti = tensorIndex(i);
    for (unsigned j = 0; j < d; ++j)
        res.row(j) =  m_bases[j]->support( ti[j] );
    return res;
}

template<unsigned d, class Basis_t >
void gsTensorBasis<d,Basis_t>::evalSingle_into(unsigned i, 
                                               const gsMatrix<T> & u, 
                                               gsMatrix<T>& result) const
{
    result.resize(1,u.cols() );
    gsMatrix<T> ev;
    gsVector<unsigned, d> ti;
    ti = tensorIndex(i);

    // Evaluate univariate basis functions and compute the product
    this->m_bases[0]->evalSingle_into(ti[0], u.row(0), result);
    for (unsigned k = 1; k < d; ++k)
    {
        this->m_bases[k]->evalSingle_into( ti[k], u.row(k), ev);
        result = result.cwiseProduct(ev);
    }
}

template<unsigned d, class Basis_t >
void gsTensorBasis<d,Basis_t>::derivSingle_into(unsigned i, 
                                                const gsMatrix<T> & u, 
                                                gsMatrix<T>& result) const
{
    gsVector<unsigned, d> ti;
    ti.noalias() = tensorIndex(i);
    gsMatrix<T> ev, dev;
    result.setOnes(d, u.cols() );

    for (unsigned k = 0; k != d; ++k)
    {
        m_bases[k]->evalSingle_into ( ti[k], u.row(k), ev  );
        m_bases[k]->derivSingle_into( ti[k], u.row(k), dev );
 
        result.row(k)            = result.row(k).cwiseProduct(dev);
        result.topRows(k)        = result.topRows(k) * ev.asDiagonal();
        // for (unsigned j = 0; j != k; ++j)
        //     result.row(j) = result.row(j).cwiseProduct(ev);
        result.bottomRows(d-k-1) = result.bottomRows(d-k-1) * ev.asDiagonal();
        // for (unsigned j = k+1; j != d; ++j)
        //     result.row(j) = result.row(j).cwiseProduct(ev);
    }
}

template<unsigned d, class Basis_t >
void gsTensorBasis<d,Basis_t>::deriv2Single_into(unsigned i, 
                                                 const gsMatrix<T> & u, 
                                                 gsMatrix<T>& result) const
{
    gsVector<unsigned, d> ti;
    ti.noalias() = tensorIndex(i);
    gsMatrix<T> ev[d], dev[d], ddev;
    result.setOnes( d*(d + 1)/2, u.cols() );

    // Precompute values and first derivatives
    for (unsigned k = 0; k != d; ++k)
    {
        m_bases[k]->evalSingle_into ( ti[k], u.row(k),  ev[k]  );
        m_bases[k]->derivSingle_into( ti[k], u.row(k), dev[k]  );
    }

    int c = d;
    for (unsigned k = 0; k != d; ++k)
    {
        // Pure second derivatives
        m_bases[k]->deriv2Single_into( ti[k], u.row(k), ddev );
        result.row(k)                = result.row(k).cwiseProduct(ddev);

        // Multiply values to all other second derivatives
        result.topRows(k)            = result.topRows(k) * ev[k].asDiagonal();
        result.middleRows(k+1,d-k-1) = result.middleRows(k+1, d-k-1) * ev[k].asDiagonal();

        // Second mixed derivatives
        for (unsigned l = k+1; l != d; ++l)
        {
            // Multiply with k-th first derivative
            result.row(c)     = result.row(c).cwiseProduct(dev[k]);
            // Multiply with l-th first derivative
            result.row(c)     = result.row(c).cwiseProduct(dev[l]);

            // Multiply with values
            for (unsigned r = 0; r != k; ++r)
                result.row(c) = result.row(c).cwiseProduct(ev[r]);
            for (unsigned r = k+1; r != l; ++r)
                result.row(c) = result.row(c).cwiseProduct(ev[r]);
            for (unsigned r = l+1; r != d; ++r)
                result.row(c) = result.row(c).cwiseProduct(ev[r]);
            c++;
        }
    }
}

template<unsigned d, class Basis_t >
void gsTensorBasis<d,Basis_t>::eval_into(const gsMatrix<T> & u, 
                                         gsMatrix<T>& result) const
{
    GISMO_ASSERT( u.rows() == d, 
                  "Attempted to evaluate the tensor-basis on points with the wrong dimension" );

    /// \todo more efficient ?

    gsMatrix<T> ev[d];
    gsVector<unsigned, d> v, size;

    // Evaluate univariate basis functions
    unsigned nb = 1;
    for (unsigned i = 0; i < d; ++i)
    {
        m_bases[i]->eval_into( u.row(i), ev[i] );
        nb *= ev[i].rows();
        size[i] = ev[i].rows();
    }
  
    // initialize result
    result.resize( nb, u.cols() );

    // iterate over all tensor product basis functions
    v.setZero();
    unsigned r = 0;
    do {
        // Multiply BSplines to get the value of basis function v
        result.row( r )=  ev[0].row( v(0) );
        for ( unsigned i=1; i<d; ++i)
            result.row( r )= result.row( r ).cwiseProduct( ev[i].row( v(i) ) ) ;
      
        ++r ;
    } while (nextLexicographic(v, size));
};

template<unsigned d, class Basis_t >
void gsTensorBasis<d,Basis_t>::eval_into(const gsMatrix<T> & u, 
                                         const gsMatrix<T> & coefs, 
                                         gsMatrix<T>& result ) const
{

#if (0)
    //~ version using the tensor structure; uses more memory, check
    result.resize( coefs.cols(), u.cols() ) ;

    unsigned sz = this->size()/m_bases[0]->size();
    unsigned n  = coefs.cols();
    gsAsMatrix<T> cc(coefs.data(), m_bases[0]->size(),  n * sz);
    gsMatrix<T> e0,e1; 

    for ( index_t j=0; j< u.cols() ; j++ ) // for all points (columns of u)
    { // to try:: all points at once instead
        gsMatrix<T> uu = u.col(j);
        //gsDebug<< "uu : \n"<< uu <<std::endl;
        //gsDebug<< "coefs: \n"<< cc <<std::endl;
        
        m_bases[0]->eval_into(uu.row(0), cc, e0);
        
        for ( unsigned i=1; i< d ; ++i )
        {
            sz /= m_bases[i]->size();
            e0.resize( m_bases[i]->size(),  n * sz );
            //gsDebug<< "e0 : \n"<< e0  <<"\n";
            m_bases[i]->eval(uu.row(i), e0, e1);
            //gsDebug<< "e1 : \n"<< e1 <<"\n";
            std::swap(e0, e1);
        }
        result->col(j) = e0;
    }
    //gsDebug<< "final: \n"<< *result  <<std::endl;
#endif
  
    // Default version - linear combination of basis functions
    result.resize( coefs.cols(), u.cols() ) ;
    gsMatrix<T> B ;
    gsMatrix<unsigned> ind;
    
    // "eval" of gsTensorBasis
    gsTensorBasis::eval_into(u, B);   // col j = nonzero basis functions at column point u(..,j)
    gsTensorBasis::active_into(u,ind);// col j = indices of active functions at column point u(..,j)
    
    for ( index_t j=0; j< u.cols() ; j++ ) // for all points (columns of u)
    {
        result.col(j) = coefs.row( ind(0,j) ) * B(0,j);
        for ( index_t i = 1; i < ind.rows(); ++i ) // for all non-zero basis functions
            result.col(j) += coefs.row( ind(i,j) ) * B(i,j); 
    }
}


template<unsigned d, class Basis_t >
void gsTensorBasis<d,Basis_t>::deriv_into(const gsMatrix<T> & u, 
                                          gsMatrix<T>& result) const
{
    gsMatrix<T> values[d];

    gsVector<unsigned, d> v, size;

    unsigned nb = 1;
    for (unsigned i = 0; i < d; ++i)
    {
        m_bases[i]->evalAllDers_into( u.row(i), 1, values[i]); // evaluate basis functions and their first derivatives

        const int num_i = values[i].rows() / 2;            // each basis function has a value and a derivative
        nb *= num_i;
        size[i] = num_i;
    }

    result.resize( d*nb, u.cols() );

    v.setZero();
    unsigned r = 0;
    do {
        for ( unsigned k=0; k<d; ++k)
        {
            const index_t rownum = r*d + k;
            result.row(rownum)  =  values[k].row( size[k] + v(k) );                       // derivative w.r.t. k-th variable
            for ( unsigned i=0; i<k; ++i)
                result.row(rownum).array() *= values[i].row( v(i) ).array();
            for ( unsigned i=k+1; i<d; ++i)
                result.row(rownum).array() *= values[i].row( v(i) ).array();
        }
        ++r ;
    } while (nextLexicographic(v, size));
}


template<unsigned d, class Basis_t >
void gsTensorBasis<d,Basis_t>::evalAllDers_into(const gsMatrix<T> & u, int n, 
                                                gsMatrix<T>& result ) const
{
    GISMO_ASSERT( n<=2, "gsTensorBasis::evalAllDers() not implemented for n > 1." );

    gsMatrix<T> values[d];

    gsVector<unsigned, d> v, size;

    unsigned nb = 1;
    for (unsigned i = 0; i < d; ++i)
    {
        // evaluate basis functions/derivatives
        m_bases[i]->evalAllDers_into( u.row(i), n, values[i] ); 
      
        // each basis function has n+1 values
        const int num_i = values[i].rows() / (n+1);
        nb *= num_i;
        size[i] = num_i;
    }
    if (n==2)
        result.resize((1 + 2*d + (d*(d-1))/2)*nb, u.cols());
    else
        result.resize((1 + n*d)*nb, u.cols());
    // iterate over all tensor product basis functions
    v.setZero();
    unsigned r = 0;
    do {
        // Multiply BSplines to get the value of basis function v
        result.row( r )=  values[0].row( v(0) );
        for ( unsigned i=1; i<d; ++i)
            result.row( r ).array() *= values[i].row( v(i) ).array();
    
        ++r ;
    } while (nextLexicographic(v, size));
  
    // iterate again and write derivatives
    if ( n>=1)
    {
        v.setZero();
        do {
            for ( unsigned k=0; k<d; ++k)
            {
                const index_t rownum = nb + (r-nb)*d + k;
                // derivative w.r.t. k-th variable
                result.row(rownum)  =  values[k].row( size[k] + v(k) );
                for ( unsigned i=0; i<k; ++i)
                    result.row(rownum).array() *= values[i].row( v(i) ).array();
                for ( unsigned i=k+1; i<d; ++i)
                    result.row(rownum).array() *= values[i].row( v(i) ).array();
            }
            ++r ;
        } while (nextLexicographic(v, size));
    }
    if (n==2)
    {
        gsMatrix<T> result_2der;
        deriv2_into(u,result_2der);
        result.block((1 + d)*nb,0,((d+(d*(d-1))/2)*nb),u.cols()).noalias() = result_2der;
    }
}


template<unsigned d, class Basis_t >
void gsTensorBasis<d,Basis_t>::deriv_into(const gsMatrix<T> & u, 
                                          const gsMatrix<T> & coefs, 
                                          gsMatrix<T>& result ) const 
{
    gsBasis<T>::deriv_into(u, coefs, result);
}

template<unsigned d, class Basis_t >
void gsTensorBasis<d,Basis_t>::deriv2_into(const gsMatrix<T> & u, 
                                           gsMatrix<T>& result ) const 
{
    gsMatrix<T> ev[d];
    gsMatrix<T> d1[d];
    gsMatrix<T> d2[d];

    gsVector<unsigned, d> v, size;

    unsigned nb = 1;
    for (unsigned i = 0; i < d; ++i)
    {
        // TODO: write this in terms of evalAllDers()
        m_bases[i]->eval_into   ( u.row(i), ev[i] );
        m_bases[i]->deriv_into  ( u.row(i), d1[i] );
        m_bases[i]->deriv2_into ( u.row(i), d2[i] );

        nb *= ev[i].rows();
        size[i] = ev[i].rows();
    }

    const unsigned stride = d + d*(d-1)/2;

    result.resize( stride*nb, u.cols() );

    v.setZero();
    unsigned r = 0; // r is a local index of a basis function times the stride
    do {
        unsigned m=d;
        for ( unsigned k=0; k<d; ++k)// First compute the pure second derivatives
        {
            result.row( r + k )  =  d2[k].row( v(k) ) ;// pure 2nd derivate w.r.t. k-th variable
            for ( unsigned i=0; i<k; ++i)
                result.row( r + k )= result.row( r + k ).cwiseProduct( ev[i].row( v(i) ) ) ;
            for ( unsigned i=k+1; i<d; ++i)
                result.row( r + k )= result.row( r + k ).cwiseProduct( ev[i].row( v(i) ) ) ;
          
            for ( unsigned l=k+1; l<d; ++l) // Then all mixed derivatives follow in lex order
            {
                result.row( r + m ) =
                    d1[k].row( v(k) ).cwiseProduct( d1[l].row( v(l) ) );
                for ( unsigned i=0; i<k; ++i)
                    result.row( r + m )= 
                        result.row( r + m ).cwiseProduct( ev[i].row( v(i) ) ) ;
                for ( unsigned i=k+1; i<l; ++i)
                    result.row( r + m )= 
                        result.row( r + m ).cwiseProduct( ev[i].row( v(i) ) ) ;
                for ( unsigned i=l+1; i<d; ++i)
                    result.row( r + m )= 
                        result.row(r + m).cwiseProduct( ev[i].row( v(i) ) ) ;
                ++m;
            }
        }
    
        r+= stride;
    } while (nextLexicographic(v, size));
}

template <unsigned d, class Basis_t>
void gsTensorBasis<d,Basis_t>::uniformRefine_withCoefs(gsMatrix<T>& coefs, int numKnots)
{
    // Simple implementation: get the transfer matrix and apply it.
    // Could be done more efficiently if needed.
    gsSparseMatrix<T, RowMajor> transfer;
    this->uniformRefine_withTransfer(transfer, numKnots);
    coefs = transfer * coefs;
}


template <unsigned d, class Basis_t>
void gsTensorBasis<d,Basis_t>::uniformRefine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, int numKnots)
{
    gsSparseMatrix<T,RowMajor> B[d];

    // refine component bases and obtain their transfer matrices
    for (unsigned i = 0; i < d; ++i)
    {
        m_bases[i]->uniformRefine_withTransfer( B[i], numKnots );
    }

    tensorCombineTransferMatrices<d, T>( B, transfer );
}

/*
 * //Note: MSVC won't resolve this if defined outside the class
template<unsigned d, class Basis_t >
typename gsBasis<typename gsTensorBasis<d,Basis_t>::T>::domainIter
gsTensorBasis<d,Basis_t>::makeDomainIterator() const
{
    return typename gsBasis<T>::domainIter(new gsTensorDomainIterator<T, d>(*this));
}


template<class T>
template<unsigned d, class Basis_t >
std::auto_ptr<gsDomainIterator<T> >
//typename memory::auto_ptr<gsDomainIterator<typename gsTensorBasis<d,Basis_t>::T> >
//typename gsBasis<typename gsTensorBasis<d,Basis_t>::T>::domainIter
gsTensorBasis<d,Basis_t>::makeDomainIterator(const boxSide & s) const
{
    return typename gsBasis<T>::domainIter(new gsTensorDomainBoundaryIterator<T, d>(*this,s));
}
*/


template <unsigned d, class Basis_t>
gsGeometry<typename gsTensorBasis<d,Basis_t>::T> * 
gsTensorBasis<d,Basis_t>::interpolate(gsMatrix<T> const& vals) const
{
    std::vector<gsMatrix<T> > grid(d);

    for (unsigned i = 0; i < d; ++i) // for all coordinate bases
        m_bases[i]->anchors_into(grid[i]);
    
    return interpolateGrid(vals,grid);
}


template <unsigned d, class Basis_t>
gsGeometry<typename gsTensorBasis<d,Basis_t>::T> * 
gsTensorBasis<d,Basis_t>::interpolateGrid(gsMatrix<T> const& vals,
                                          std::vector<gsMatrix<T> >const& grid) const
{
    GISMO_ASSERT (this->size() == vals.cols(), 
                  "Expecting as many values as the number of basis functions." );

    const index_t n  = vals.rows();
    const int sz = this->size();
    int sz_i, r_i;

    // Note: algorithm relies on col-major matrices    
    gsMatrix<T> q0, q1;

    //Note: Sparse LU might fail for rank deficient Cmat
    Eigen::SparseLU<gsSparseMatrix<T>, Eigen::COLAMDOrdering<index_t> >  solver;
    gsSparseMatrix<T> Cmat;

    q0 = vals.transpose();

    for (unsigned i = 0; i < d; ++i) // for all coordinate bases
    {
        // Re-order right-hand sides
        sz_i = m_bases[i]->size();
        r_i  = sz / sz_i;
        q0.resize(sz_i, n * r_i);

        // Solve for i-th coordinate basis
        m_bases[i]->collocationMatrix(grid[i], Cmat);
        solver.compute(Cmat); 
        if ( solver.info() != Eigen::Success )
        {
            gsWarn<< "Failed LU decomposition for:\n";//<< Cmat.toDense() <<"\n";
            gsWarn<< "Points:\n"<< grid[i] <<"\n";
            gsWarn<< "Knots:\n"<< m_bases[i]->detail() <<"\n";
            return 0;
        }

        // Transpose solution component-wise
        q1.resize(r_i, n * sz_i);
        for ( index_t k = 0; k!=n; ++k)
            q1.middleCols(k*sz_i, sz_i) = solver.solve(q0.middleCols(k*r_i,r_i)).transpose();

        q1.swap( q0 ); // move solution as next right-hand side
    }

    q0.resize(sz, n);
    return this->makeGeometry( give(q0) );
}


//template <unsigned d, class Basis_t>
//gsDomain<typename gsTensorBasis<d,Basis_t>::T> * gsTensorBasis<d,Basis_t>::makeDomain() const 
//{
//  return new gsTensorDomain<T>();
//} 
 

} // namespace gismo
