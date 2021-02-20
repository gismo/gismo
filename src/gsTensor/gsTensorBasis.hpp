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
#include <gsCore/gsGeometry.h>
//#include <gsUtils/gsSortedVector.h>


#if defined __INTEL_COMPILER
#pragma warning (disable : 175 ) /* disable subscript out of range warning */
#endif

namespace gismo
{

template<short_t d, class T>
gsTensorBasis<d,T>::gsTensorBasis( Basis_t* x,  Basis_t* y)
{
    GISMO_ASSERT( d==2, "gsTensorBasis error: wrong dimension." );

    if ( x->dim() == 1 && y->dim() == 1 )
    {
        m_bases[0] = x;
        m_bases[1] = y;
    } else 
        GISMO_ERROR("gsTensorBasis error: Spaces must be of topological dimension 1.");
}


template<short_t d, class T>
gsTensorBasis<d,T>::gsTensorBasis(Basis_t* x, Basis_t* y, Basis_t* z)
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

template<short_t d, class T>
gsTensorBasis<d,T>::gsTensorBasis(Basis_t* x, Basis_t* y, Basis_t* z, Basis_t* w)
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


/*
template<short_t d, class T>
gsTensorBasis<d,T>::gsTensorBasis( std::vector<Basis_t*> const & bb )
{ 
    GISMO_ASSERT( d==bb.size(), "gsTensorBasis error: wrong number of bases given ("
                  << bb.size()<< ", expected "<< d );

    for (unsigned i = 0; i < d; ++i)
        m_bases[i] = bb[i];
}
*/

/// Copy Constructor
template<short_t d, class T>
gsTensorBasis<d,T>::gsTensorBasis( const gsTensorBasis & o)
: gsBasis<T>(o)
{
    for (short_t i = 0; i < d; ++i)
        m_bases[i] = o.m_bases[i]->clone().release();
}


template<short_t d, class T>
gsTensorBasis<d,T>& gsTensorBasis<d,T>::operator=( const gsTensorBasis & o)
{
    if ( this == &o )
        return *this;
    gsBasis<T>::operator=(o);

    for (short_t i = 0; i < d; ++i)
    {
        delete m_bases[i];
        m_bases[i] = o.m_bases[i]->clone().release();
    }
    return *this;
}


template<short_t d, class T>
void gsTensorBasis<d,T>::anchors_into(gsMatrix<T>& result) const
{
    gsMatrix<T> gr[d];
    gsVector<unsigned, d> v, size;
    result.resize( d, this->size() );

    for (short_t i = 0; i < d; ++i)
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

template<short_t d, class T>
void gsTensorBasis<d,T>::anchor_into(index_t i, gsMatrix<T>& result) const
{
    gsVector<index_t, d> ti = tensorIndex(i);

    gsMatrix<T> gr;
    result.resize(d, 1);

    for (short_t l = 0; l < d; ++l)
    {
        m_bases[l]->anchor_into(ti[l], gr);
        result(l,0) = gr.value(); 
    }
}

template<short_t d, class T>
void gsTensorBasis<d,T>::connectivity(const gsMatrix<T> & nodes,
                                            gsMesh<T> & mesh) const
{
    const index_t sz  = size();
    GISMO_ASSERT( nodes.rows() == sz, "Invalid input.");

    // Add vertices
    for(index_t i = 0; i< sz; ++i )
        mesh.addVertex( nodes.row(i).transpose() );

    // Starting from lower corner
    const gsVector<index_t,d> low = gsVector<index_t,d>::Zero();

    // Last tensor-index
    gsVector<index_t, d> end;
    for (short_t i = 0; i < d; ++i)
        end(i) = this->size(i) - 1;

    unsigned k, s;
    gsVector<index_t,d> v, upp;
    for (short_t i = 0; i < d; ++i) // For all axes
    {
        s      = stride(i);
        v      = low;
        upp    = end;
        upp[i] = 0; // suppress to face v[i]==0

        do // Insert all edges normal to axis i 
        {
            k = index(v);
            for (index_t j = 0; j != end[i]; ++j)
            {
                mesh.addEdge(k, k + s );
                k += s;
            }
        }
        while ( nextCubePoint(v, low, upp) );
    }
}

template<short_t d, class T>
void gsTensorBasis<d,T>::active_into(const gsMatrix<T> & u, gsMatrix<index_t>& result) const
{
    //gsWarn<<"genericActive "<< *this;

    gsMatrix<index_t> act[d];
    gsVector<index_t, d> v, size;
 
    // Get component active basis functions
    index_t nb = 1;
    for (short_t i = 0; i < d; ++i)
    {
        m_bases[i]->active_into(u.row(i), act[i]);
        size[i] = act[i].rows();
        nb     *= size[i];
    }

    result.resize( nb, u.cols() );
   
    // iterate over all tensor product active functions
    unsigned r = 0;
    v.setZero();
    do {
        // Fill with active bases indices
        for ( index_t j=0; j<u.cols(); ++j)
        {
            nb = act[d-1]( v(d-1) ,j) ;//compute global index in the tensor product
            for ( short_t i=d-2; i>=0; --i ) // to do: strides
                nb = nb * m_bases[i]->size() + act[i](v(i), j) ;
            result(r, j) = nb;
        }
        ++r ;
    } while (nextLexicographic(v, size));
}

template<short_t d, class T>
bool gsTensorBasis<d,T>::isActive(const index_t i, const gsVector<T>& u) const
{
    GISMO_ASSERT( u.rows() == static_cast<index_t>(d), "Invalid input.");
    const gsVector<index_t, d> ti = tensorIndex(i);
    for (short_t k = 0; k < d; ++k)
        if (  ! m_bases[k]->isActive(ti[k], u.row(k)) )
            return false;
    return true;
}

template<short_t d, class T>
gsMatrix<index_t> gsTensorBasis<d,T>::coefSlice(short_t dir, index_t k) const
{
    GISMO_ASSERT( dir>=0 &&  dir < this->dim(), "Invalid slice direction requested" );
    GISMO_ASSERT( k >=0 &&  k < this->size(dir), "Invalid slice position requested" );

    index_t sliceSize = 1, r = 0;
    gsVector<index_t, d> low, upp;

    // Set low and upp to point to the range of indices
    for (short_t i = 0; i < d; ++i)
    {
        sliceSize *= this->size(i); // or trueSize(i) ?
        low(r) = 0;
        upp(r++) = this->size(i); // or trueSize(i) ?
    }
    sliceSize /= upp(dir);
    upp(dir) = k + 1;
    low(dir) = k;

    gsMatrix<index_t> res(sliceSize,1);

    // iterate over all tensor product basis indices
    gsVector<index_t,d> v = low;
    r = 0;
    do {
        res(r++,0) = this->index(v);
    } while ( nextLexicographic(v, low, upp) );

    return res;
}


template<short_t d, class T>
gsMatrix<index_t> gsTensorBasis<d,T>::allBoundary() const
{
    gsMatrix<index_t> bd;
    std::set<index_t> bdofs;

    for (short_t k = 0; k != d; ++k)
    {
        bd = this->coefSlice(k, 0);
        for (index_t i = 0; i < bd.size(); ++i)
            bdofs.insert( bd(i) );

        bd = this->coefSlice(k, size(k) - 1);
        for (index_t i = 0; i < bd.size(); ++i)
            bdofs.insert( bd(i) );
    }

    return makeMatrix<index_t>(bdofs.begin(), static_cast<index_t>(bdofs.size()), 1 );
}


template<short_t d, class T>
gsMatrix<index_t> gsTensorBasis<d,T>::boundaryOffset(boxSide const& s,index_t offset) const
{
    //get m_bases index and start or end case
    short_t k = s.direction();
    bool r = s.parameter();
    GISMO_ASSERT(offset < size(k),"Offset cannot be bigger than the amount of basis functions orthogonal to Boxside s!");
    return (this->coefSlice(k, (r ? size(k) - 1 -offset : offset) ));
}

template<short_t d, class T>
index_t gsTensorBasis<d,T>::functionAtCorner(boxCorner const & c) const
{
    gsVector<bool> position(d);
    c.parameters_into(d, position);
    
    index_t index = 0;
    index_t str   = 1;
    
    for(short_t i = 0; i!=d; ++i)
    {
        const index_t sz_i = size(i);
        if ( position[i] )
            index+= str * ( sz_i - 1 );
        str *= sz_i;
    }

    return index;
}

/*
template<short_t d, class T>
void gsTensorBasis<d,T>::boundary_into(boxSide const & s, gsMatrix<int> & bstruct, gsMatrix<unsigned>& result) const
{
    //get m_bases index and start or end case
    index_t k = s.direction();
    index_t r = s.parameter();

    // Compute the structure of the boundary dofs
    bstruct.resize(d-1);
    index_t c = 0;
    for (index_t k = 0; k<d; ++k )
    {
        if ( k == s ) continue;
        bSize[c] = m_bases[k]->size();
        c++;
    }

    return this->slice(k, (r ? size(k) - 1 : 0) ).release();
}
//*/

/*
template <short_t d, class BB, class B>
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
//*/

template<short_t d, class T>
void
gsTensorBasis<d,T>::getComponentsForSide(boxSide const& s, std::vector<Basis_t*> & rr) const
{   
    index_t dir = s.direction( );

    rr.clear();
    rr.reserve( d - 1 );
    for ( short_t i=0; i < d; ++i )
        if (i != dir)
            rr.push_back(m_bases[i]->clone().release());
}


template<short_t d, class T>
gsMatrix<T> gsTensorBasis<d,T>::support() const
{
    gsMatrix<T> res(d,2);
    for (short_t i = 0; i < d; ++i)
        res.row(i) =  m_bases[i]->support();
    return res;
}

template<short_t d, class T>
gsMatrix<T>
gsTensorBasis<d,T>::support(const index_t & i) const
{
    gsMatrix<T> res(d,2);
    gsVector<index_t, d> ti = tensorIndex(i);
    for (short_t j = 0; j < d; ++j)
        res.row(j) =  m_bases[j]->support( ti[j] );
    return res;
}

template<short_t d, class T>
void gsTensorBasis<d,T>::evalSingle_into(index_t i,
                                               const gsMatrix<T> & u,
                                               gsMatrix<T>& result) const
{
    result.resize(1,u.cols() );
    gsMatrix<T> ev;
    gsVector<index_t, d> ti;
    ti = tensorIndex(i);

    // Evaluate univariate basis functions and compute the product
    this->m_bases[0]->evalSingle_into(ti[0], u.row(0), result);
    for (short_t k = 1; k < d; ++k)
    {
        this->m_bases[k]->evalSingle_into( ti[k], u.row(k), ev);
        result = result.cwiseProduct(ev);
    }
}

template<short_t d, class T>
void gsTensorBasis<d,T>::derivSingle_into(index_t i,
                               const gsMatrix<T> & u,
                                      gsMatrix<T>& result) const
{
    gsVector<index_t, d> ti;
    ti.noalias() = tensorIndex(i);
    gsMatrix<T> ev, dev;
    result.setOnes(d, u.cols() );

    for (short_t k = 0; k != d; ++k)
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

template<short_t d, class T>
void gsTensorBasis<d,T>::deriv2Single_into(index_t i,
                                const gsMatrix<T> & u,
                                       gsMatrix<T>& result) const
{
    gsVector<index_t, d> ti;
    ti.noalias() = tensorIndex(i);
    gsMatrix<T> ev[d], dev[d], ddev;
    result.setOnes( d*(d + 1)/2, u.cols() );

    // Precompute values and first derivatives
    for (short_t k = 0; k != d; ++k)
    {
        m_bases[k]->evalSingle_into ( ti[k], u.row(k),  ev[k]  );
        m_bases[k]->derivSingle_into( ti[k], u.row(k), dev[k]  );
    }

    index_t c = d;
    for (short_t k = 0; k != d; ++k)
    {
        // Pure second derivatives
        m_bases[k]->deriv2Single_into( ti[k], u.row(k), ddev );
        result.row(k)                = result.row(k).cwiseProduct(ddev);

        // Multiply values to all other second derivatives
        result.topRows(k)            = result.topRows(k) * ev[k].asDiagonal();
        result.middleRows(k+1,d-k-1) = result.middleRows(k+1, d-k-1) * ev[k].asDiagonal();

        // Second mixed derivatives
        for (short_t l = k+1; l != d; ++l)
        {
            // Multiply with k-th first derivative
            result.row(c)     = result.row(c).cwiseProduct(dev[k]);
            // Multiply with l-th first derivative
            result.row(c)     = result.row(c).cwiseProduct(dev[l]);

            // Multiply with values
            for (short_t r = 0; r != k; ++r)
                result.row(c) = result.row(c).cwiseProduct(ev[r]);
            for (short_t r = k+1; r != l; ++r)
                result.row(c) = result.row(c).cwiseProduct(ev[r]);
            for (short_t r = l+1; r != d; ++r)
                result.row(c) = result.row(c).cwiseProduct(ev[r]);
            c++;
        }
    }
}

template<short_t d, class T>
void gsTensorBasis<d,T>::eval_into(const gsMatrix<T> & u,
                                         gsMatrix<T>& result) const
{
    GISMO_ASSERT( u.rows() == d, 
                  "Attempted to evaluate the tensor-basis on points with the wrong dimension" );

    /// \todo more efficient ?

    gsMatrix<T> ev[d];
    gsVector<unsigned, d> v, size;

    // Evaluate univariate basis functions
    unsigned nb = 1;
    for (short_t i = 0; i < d; ++i)
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
        for ( short_t i=1; i<d; ++i)
            result.row( r )= result.row( r ).cwiseProduct( ev[i].row( v(i) ) ) ;
      
        ++r ;
    } while (nextLexicographic(v, size));
};

template<short_t d, class T>
void gsTensorBasis<d,T>::eval_into(const gsMatrix<T> & u,
                                   const gsMatrix<T> & coefs,
                                         gsMatrix<T> & result ) const
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
        
        for ( short_t i=1; i< d ; ++i )
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
    // Basis_t:: eval_into(u,coefs,result);
    result.resize( coefs.cols(), u.cols() ) ;
    gsMatrix<T> B ;
    gsMatrix<index_t> ind;
    
    this->eval_into(u, B);   // col j = nonzero basis functions at column point u(..,j)
    this->active_into(u,ind);// col j = indices of active functions at column point u(..,j)
    
    for ( index_t j=0; j< u.cols() ; j++ ) // for all points (columns of u)
    {
        result.col(j) = coefs.row( ind(0,j) ) * B(0,j);
        for ( index_t i = 1; i < ind.rows(); ++i ) // for all non-zero basis functions
            result.col(j) += coefs.row( ind(i,j) ) * B(i,j); 
    }
}


template<short_t d, class T>
void gsTensorBasis<d,T>::deriv_into(const gsMatrix<T> & u,
                                          gsMatrix<T>& result) const
{
    std::vector<gsMatrix<T> > values[d];

    gsVector<unsigned, d> v, size;

    index_t nb = 1;
    for (short_t i = 0; i < d; ++i)
    {
        // evaluate basis functions and their first derivatives
        m_bases[i]->evalAllDers_into( u.row(i), 1, values[i]); 

        // number of basis functions
        const index_t num_i = values[i].front().rows();
        nb *= num_i;
        size[i] = num_i;
    }

    result.resize( d*nb, u.cols() );

    v.setZero();
    unsigned r = 0;
    do {
        for ( short_t k=0; k<d; ++k)
        {
            // derivative w.r.t. k-th variable
            const index_t rownum = r*d + k;
            result.row(rownum)  =  values[k][1].row( v(k) );
            for ( short_t i=0; i<k; ++i)
                result.row(rownum).array() *= values[i][0].row( v(i) ).array();
            for ( short_t i=k+1; i<d; ++i)
                result.row(rownum).array() *= values[i][0].row( v(i) ).array();
        }
        ++r ;
    } while (nextLexicographic(v, size));
}


template<short_t d, class T>
void gsTensorBasis<d,T>::evalAllDers_into(const gsMatrix<T> & u, int n,
                                          std::vector<gsMatrix<T> >& result) const
{
    GISMO_ASSERT(n>-2, "gsTensorBasis::evalAllDers() is implemented only for -2<n<=2: -1 means no value, 0 values only, ... " );
    if (n==-1)
    {
        result.resize(0);
        return;
    }

    std::vector< gsMatrix<T> >values[d];
    gsVector<unsigned, d> v, nb_cwise;
    result.resize(n+1);

    unsigned nb = 1;
    for (short_t i = 0; i < d; ++i)
    {
        // evaluate basis functions/derivatives
        m_bases[i]->evalAllDers_into( u.row(i), n, values[i] ); 
      
        // number of basis functions
        const index_t num_i = values[i].front().rows();
        nb_cwise[i] = num_i;
        nb         *= num_i;
    }
    
    // iterate over all tensor product basis functions
    v.setZero();
    gsMatrix<T> & vals = result[0];
    vals.resize(nb, u.cols());
    unsigned r = 0;
    do // for all basis functions
    {
        // Multiply basis functions to get the value of basis function v
        vals.row( r )=  values[0].front().row( v[0] );
        for ( short_t i=1; i!=d; ++i) // for all variables
            vals.row(r).array() *= values[i].front().row( v[i] ).array();
    
        ++r;
    } while (nextLexicographic(v, nb_cwise));
  
    // iterate again and write derivatives
    if ( n>=1)
    {
        gsMatrix<T> & der = result[1];
        der.resize(d*nb, u.cols());;
        v.setZero();
        r = 0;
        do // for all basis functions
        {
            for ( short_t k=0; k<d; ++k) // for partial derivatives
            {
                // derivative w.r.t. k-th variable // TODO Check this!
                der.row(r) = values[k][1].row( v(k) );
                for ( short_t i=0; i<k; ++i)
                    der.row(r).array() *= values[i][0].row( v(i) ).array();
                for ( short_t i=k+1; i<d; ++i)
                    der.row(r).array() *= values[i][0].row( v(i) ).array();
            ++r;
            }

        } while (nextLexicographic(v, nb_cwise));
    }

    if (n>1)
    {
        deriv2_tp( values, nb_cwise, result[2] );

        gsVector<unsigned, d> cc;
        for (int i = 3; i <=n; ++i) // for all orders of derivation
        {
            gsMatrix<T> & der = result[i];
            der.resize( nb*numCompositions(i,d), u.cols());
            v.setZero();
            
            r = 0;
            do // for all basis functions
            {
                firstComposition(i, d, cc);
                do // for all partial derivatives of order \a i
                {
                    // cc[k]: order of derivation w.r.t. variable \a k
                    der.row(r) = values[0][cc[0]].row(v[0]);
                    for (short_t k = 1; k!=d; ++k) // for all variables
                        der.row(r).array() *= values[k][cc[k]].row(v[k]).array();
                ++r;
                } while (nextComposition(cc));
            } while (nextLexicographic(v, nb_cwise));
        }

        // n==2: bubble up // <- not needed due to call of deriv2_tp
    }

}

template<short_t d, class T>
void gsTensorBasis<d,T>::deriv2_into(const gsMatrix<T> & u,
                                           gsMatrix<T> & result ) const
{
    std::vector< gsMatrix<T> >values[d];
    gsVector<unsigned, d> v, nb_cwise;

    unsigned nb = 1;
    for (short_t i = 0; i < d; ++i)
    {
        m_bases[i]->evalAllDers_into( u.row(i), 2, values[i]); 
        const int num_i = values[i].front().rows();
        nb_cwise[i] = num_i;
        nb     *= num_i;
    }

    deriv2_tp(values, nb_cwise, result);
}



template<short_t d, class T>
void gsTensorBasis<d,T>::deriv2_tp(const std::vector< gsMatrix<T> > values[],
                                   const gsVector<unsigned, d> & nb_cwise,
                                   gsMatrix<T>& result)
{
    const unsigned nb = nb_cwise.prod();
    const unsigned stride = d + d*(d-1)/2;

    result.resize( stride*nb, values[0][0].cols() );

    gsVector<unsigned, d> v;
    v.setZero();
    unsigned r = 0; // r is a local index of a basis function times the stride
    do
    {
        unsigned m = d;
        for ( short_t k=0; k<d; ++k)// First compute the pure second derivatives
        {
            index_t cur = r + k;
            result.row(cur) = values[k][2].row( v.at(k) ) ;// pure 2nd derivate w.r.t. k-th var
            for ( short_t i=0; i<k; ++i)
                result.row(cur).array() *= values[i][0].row( v.at(i) ).array();
            for ( short_t i=k+1; i<d; ++i)
                result.row(cur).array() *= values[i][0].row( v.at(i) ).array();

            for ( short_t l=k+1; l<d; ++l) // Then all mixed derivatives follow in lex order
            {
                cur = r + m;
                result.row(cur).noalias() =
                    values[k][1].row( v.at(k) ).cwiseProduct( values[l][1].row( v.at(l) ) );
                for ( short_t i=0; i<k; ++i)
                    result.row(cur).array() *= values[i][0].row( v.at(i) ).array() ;
                for ( short_t i=k+1; i<l; ++i)
                    result.row(cur).array() *= values[i][0].row( v.at(i) ).array();
                for ( short_t i=l+1; i<d; ++i)
                    result.row(cur).array() *= values[i][0].row( v.at(i) ).array() ;
                ++m;
            }
        }

        r+= stride;
    } while (nextLexicographic(v, nb_cwise));
}


template<short_t d, class T>
void gsTensorBasis<d,T>::refineElements(std::vector<index_t> const & elements)
{
    gsSortedVector<index_t> elIndices[d];
    index_t tmp, mm;
    
    // Get coordinate wise element indices
    for ( typename  std::vector<index_t>::const_iterator
              it = elements.begin(); it != elements.end(); ++it )
    {
        mm = *it;
            for (short_t i = 0; i<d; ++i )
            {
                const index_t nEl_i = m_bases[i]->numElements();
                tmp = mm % nEl_i;
                mm = (mm - tmp) / nEl_i;
                elIndices[i].push_sorted_unique(tmp);
            }
    }
    
    // Refine in each coordinate
    // Element refinement propagates along knot-lines
    for (short_t i = 0; i<d; ++i )
    {
        m_bases[i]->refineElements(elIndices[i]);
    }
}


template<short_t d, class T>
void gsTensorBasis<d,T>::uniformRefine_withCoefs(gsMatrix<T>& coefs, int numKnots, int mul)
{
    // Simple implementation: get the transfer matrix and apply it.
    // Could be done more efficiently if needed.
    gsSparseMatrix<T, RowMajor> transfer;
    this->uniformRefine_withTransfer( transfer, numKnots, mul );
    coefs = transfer * coefs;
}


template<short_t d, class T>
void gsTensorBasis<d,T>::uniformRefine_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, int numKnots, int mul)
{
    gsSparseMatrix<T,RowMajor> B[d];

    // refine component bases and obtain their transfer matrices
    for (short_t i = 0; i < d; ++i)
    {
        m_bases[i]->uniformRefine_withTransfer( B[i], numKnots, mul );
    }

    tensorCombineTransferMatrices<d, T>( B, transfer );
}


template<short_t d, class T>
void gsTensorBasis<d,T>::uniformCoarsen_withTransfer(gsSparseMatrix<T,RowMajor> & transfer, int numKnots)
{
    gsSparseMatrix<T,RowMajor> B[d];

    // refine component bases and obtain their transfer matrices
    for (short_t i = 0; i < d; ++i)
    {
        m_bases[i]->uniformCoarsen_withTransfer( B[i], numKnots );
    }

    tensorCombineTransferMatrices<d, T>( B, transfer );
}

/*
 * //Note: MSVC won't resolve this if defined outside the class
template<short_t d, class T>
typename gsBasis<T>::domainIter
gsTensorBasis<d,T>::makeDomainIterator() const
{
    return typename gsBasis<T>::domainIter(new gsTensorDomainIterator<T, d>(*this));
}


template<class T>
template<short_t d, class T>
gsDomainIterator<T>::ptr
//memory::unique_ptr<gsDomainIterator<T>
//typename gsBasis<T>::domainIter
gsTensorBasis<d,T>::makeDomainIterator(const boxSide & s) const
{
    return typename gsBasis<T>::domainIter(new gsTensorDomainBoundaryIterator<T, d>(*this,s));
}
*/


template<short_t d, class T>
typename gsGeometry<T>::uPtr
gsTensorBasis<d,T>::interpolateAtAnchors(gsMatrix<T> const& vals) const
{
    std::vector<gsMatrix<T> > grid(d);

    for (short_t i = 0; i < d; ++i) // for all coordinate bases
        m_bases[i]->anchors_into(grid[i]);
    
    return interpolateGrid(vals,grid);
}


template<short_t d, class T>
typename gsGeometry<T>::uPtr
gsTensorBasis<d,T>::interpolateGrid(gsMatrix<T> const& vals,
                                          std::vector<gsMatrix<T> >const& grid) const
{
    GISMO_ASSERT (this->size() == vals.cols(), 
                  "Expecting as many values as the number of basis functions." );

    const index_t n  = vals.rows();
    const int sz = this->size();

    // Note: algorithm relies on col-major matrices    
    gsMatrix<T, Dynamic, Dynamic, ColMajor> q0, q1;

    //Note: Sparse LU might fail for rank deficient Cmat
    typename gsSparseSolver<T>::LU  solver;
    gsSparseMatrix<T> Cmat;

    // size: sz x n
    q0 = vals.transpose();

    for (short_t i = 0; i < d; ++i) // for all coordinate bases
    {
        // Re-order right-hand sides
        const index_t sz_i = m_bases[i]->size();
        const index_t r_i  = sz / sz_i;
        q0.resize(sz_i, n * r_i);

        // Solve for i-th coordinate basis
        m_bases[i]->collocationMatrix(grid[i], Cmat);
        solver.compute(Cmat); 
        #ifndef NDEBUG
        if ( solver.info() != Eigen::Success )
        {
            gsWarn<< "Failed LU decomposition for:\n";//<< Cmat.toDense() <<"\n";
            gsWarn<< "Points:\n"<< grid[i] <<"\n";
            gsWarn<< "Knots:\n"<< m_bases[i]->detail() <<"\n";
            return typename gsGeometry<T>::uPtr();
        }
        #endif
        // Transpose solution component-wise
        q1.resize(r_i, n * sz_i);
        for ( index_t k = 0; k!=n; ++k)
            q1.middleCols(k*sz_i, sz_i) = solver.solve(q0.middleCols(k*r_i,r_i)).transpose();

        q1.swap( q0 ); // move solution as next right-hand side
    }

    q0.resize(sz, n);
    return this->makeGeometry( give(q0) );
}


template<short_t d, class T>
void gsTensorBasis<d,T>::matchWith(const boundaryInterface & bi,
                                   const gsBasis<T> & other,
                                   gsMatrix<index_t> & bndThis,
                                   gsMatrix<index_t> & bndOther) const
{
    if ( const Self_t * _other = dynamic_cast<const Self_t*>(&other) )
    {
        // Grab the indices to be matched
        bndThis = this->boundary( bi.first() .side() );
        bndOther= _other->boundary( bi.second().side() );
        GISMO_ASSERT( bndThis.rows() == bndOther.rows(),
                      "Input error, sizes do not match: "
                      <<bndThis.rows()<<"!="<<bndOther.rows() );
        if (bndThis.size() == 1) return;

        // Get interface data
        const index_t s1 = bi.first() .direction();
        const index_t s2 = bi.second().direction();
        const gsVector<bool>    & dirOr = bi.dirOrientation();
        const gsVector<index_t> & bMap  = bi.dirMap();
        
        // Compute the tensor structure of bndThis
        gsVector<index_t>  bSize(d-1);
        index_t c = 0;
        for (short_t k = 0; k<d; ++k )
        {
            if ( k == s1 )
                continue;
            bSize[c] = this->component(k).size();
            c++;
        }

        // Apply flips to bndThis and bndOther so that they have the
        // same orientation
        gsVector<index_t>  bPerm(d-1);
        c = 0;
        for (short_t k = 0; k<d; ++k )
        {
            if ( k == s1 ) // skip ?
                continue;
            
            if ( ! dirOr[k] ) // flip ?
                flipTensorVector(c, bSize, bndThis);
            
            bPerm[c] = ( bMap[k] < s2 ? bMap[k] : bMap[k]-1 );
            c++;
        }
    
        // Apply permutation to bndThis and bndOther so that they
        // finally match on both sides
        permuteTensorVector<index_t,-1>(bPerm, bSize, bndThis);

        return;
    }
    
    gsWarn<<"Cannot match with "<< other <<"\n";
}

template<short_t d, class T>
T gsTensorBasis<d, T>::getMinCellLength() const
{
    T h = 0;
    for (short_t i = 0; i < d; ++i)
    {
        const T sz = m_bases[i]->getMinCellLength();
        if ( sz < h || h == 0 ) h = sz;
    }
    return h;
}

template<short_t d, class T>
T gsTensorBasis<d, T>::getMaxCellLength() const
{
    T h = 0;
    for (short_t i = 0; i < d; ++i)
    {
        const T sz = m_bases[i]->getMaxCellLength();
        if ( sz > h ) h = sz;
    }
    return h;
}

template<short_t d, class T>
gsMatrix<T> gsTensorBasis<d,T>::elementInSupportOf(index_t j) const
{
    const gsVector<index_t, d> ti = tensorIndex(j);
    gsMatrix<T> el, res(d,2);
    for (short_t i = 0; i < d; ++i)
    {
        el = m_bases[i]->elementInSupportOf(ti[i]);
        res.row(i) = el;
    }
    return res;
}
    
//template<short_t d, class T>
//gsDomain<T> * gsTensorBasis<d,T>::makeDomain() const
//{
//  return new gsTensorDomain<T>();
//} 
 

} // namespace gismo
