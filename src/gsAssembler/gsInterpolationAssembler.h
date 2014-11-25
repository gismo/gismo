 
#pragma once

#include <iostream>
#include <map>
//#include <array> // c++11
//#include <tr1/tuple>

#include <gsAssembler/gsAssembler.h>
#include <gsPde/gsPoissonPde.h>

namespace gismo
{

/** 
    Assembler class using interpolation method
*/

template<class T>
class gsInterpolationAssembler  : public gsAssembler<T>
{
public:
    typedef std::pair<int, int> indices_t;
//typedef std::tr1::tuple<int, int, int> triplet_t;
// Usage
//tuplet_t t(1, 2, 2);
// std::tr1::get<0>(t);std::tr1::get<1>(t);std::tr1::get<2>(t);
// std::tr1::get<0>(t) = 3;
public:

    /// Default empty constructor
    gsInterpolationAssembler() 
        : m_geometryFactor(NULL), m_rhsApprox(NULL), m_initialized(false)
    { }
    
    /**
     *\brief Constructs an assembler from a given geometry
     *
     * The constructor saves the pointer to the geometry.
     *
     * \code 
     *  gsInterpolationAssembler(geom).stiffness(basis)
     * \endcode
     *
     * \sa stiffness(const gsBasis<T>* B)
     */  
    gsInterpolationAssembler(const gsGeometry<T> & geom)
        : gsAssembler<T>(geom), m_initialized(false)
    { 
        m_d = geom.parDim();
        this->m_geometryFactor = NULL; 
        this->m_rhsApprox = NULL;
    }

    ~gsInterpolationAssembler() 
    { 
	if ( this->m_geometryFactor )
	    delete this->m_geometryFactor;
	if ( this->m_rhsApprox )
	    delete this->m_rhsApprox;
    }
    
public:

    void setGeometry(const gsGeometry<T> & geom)
    { 
        if (this->m_geometry != &geom )
        {
            this->m_geometry = &geom;
            m_d              = geom.parDim();
            m_initialized    = false;
        }
    } 

    virtual gsSparseSystem<T>
    assemble( const gsBasis<T>& m_basis, const gsDofMapper& mapper,
              const gsMatrix<T> & ddof, const gsPde<T> & pde, int patchIndex=0)
    {
        const gsPoissonPde<T> *poisson = dynamic_cast<const gsPoissonPde<T>*>(&pde);
        if (poisson)
            return assemblePoisson(m_basis, mapper, ddof, *poisson, patchIndex);

        GISMO_ERROR("Unknown PDE type in assemble()");
    }

    /// Assembler
    gsSparseSystem<T> assemblePoisson( const gsBasis<T>& m_basis, 
                                       const gsDofMapper& mapper, 
                                       const gsMatrix<T> & ddof, 
                                       const gsPoissonPde<T> & pde, 
                                       int patchIndex=0);

    /// Assembler  implementation with statically sized temporaries where possible
    /// Assumes an open uniform knot-vector
    template <int D>
    gsSparseSystem<T> assemblePoisson_impl( const gsBasis<T>& m_basis, 
                                            const gsDofMapper& mapper, 
                                            const gsMatrix<T> & ddof, 
                                            const gsPoissonPde<T> & pde, 
                                            int patchIndex);

////////////////////////// Individual Assemblers without initial conditions
    
    gsSparseMatrix<T> * stiffness(const gsBasis<T>& B);
    gsSparseMatrix<T> * stiffness2(const gsBasis<T>& B);
    
    gsVector<T> * moments( const gsBasis<T>& B, gsFunction<T> const & f);
    
    void computeStiffnessFactor(const gsBasis<T>& m_basis, const gsFunction<T> & rhs) const ;

    void computeStiffnessFactor(const gsBasis<T>& m_basis) const ;
    void computeExactStiffnessFactor(const gsBasis<T>& m_basis) const ;

    void computeMassFactor(const gsBasis<T>& m_basis) const ;
    void computeExactMassFactor(const gsBasis<T>& m_basis) const ;

    void computeRhsApprox(const gsBasis<T>& m_basis, const gsFunction<T> & rhs) const ;

    T interpolationError(int const & n_p = 1000) const;
    
    const gsMatrix<T> & lut(int i, int k) 
    {
        GISMO_ASSERT(i<3       , "Invalid lut requested.");
        GISMO_ASSERT(k<7       , "Invalid lut requested.");
        return m_lookup[i*7+k];
    }
    
    void makeLut();

private:

    /*
      Signature
      \a c configuration
      \a m multiplicity
      \a v indices i,j flattened
    */
    inline void signature(index_t i0, 
                          index_t i1, 
                          index_t i2, 
                          int const & dir, // direction
                          int const & r,   // GeoFactor row
                          int const & s,   // GeoFactor column
                          //int & c // lut index
                          int & v,  // pattern
                          int & c,  // configuration
                          int & m,  // multiplicity
                          T & corr  // correction coefficient
        ) const
    {
        bool q0(false),q1(r==dir),q2(s==dir) ;
        
        // initial correction coeffiecient
        corr = m_h[dir];

        // Sort indices
        if (i0>i1)
        {
            std::swap(i0,i1);
            std::swap(q0,q1);
        }
        if (i1>i2)
        {
            std::swap(i1,i2);
            std::swap(q1,q2);
        }
        if (i0>i1)
        {
            std::swap(i0,i1);
            std::swap(q0,q1);
        }
        
        // Check for right boundary
        if ( int(i2)>m_sz[dir]-m_p[dir] )
        {
            // mirrior the configuration
            std::swap(i0,i2);
            std::swap(q0,q2);
            i0 = m_sz[dir] - i0;
            i1 = m_sz[dir] - i1 - i0;// i1 gets the difference
            i2 = m_sz[dir] - i2;

            // toggle derivative sign
            if (q0) corr *= -T(1);
            if (q1) corr *= -T(1);
            if (q2) corr *= -T(1);
        }
        else
        {
            i1 = i1 - i0;// i1 gets the difference
        }

        // Set pattern
        v=0;
        if (q0) 
        {
            v|=1<<2;
            corr /= m_h[dir];
        }
        if (q1)
        {
            v|=1<<1;
            corr /= m_h[dir];
        }
        if (q2)
        {
            v|=1;
            corr /= m_h[dir];
        }
        
        // Set configuration
        c = (m_p[dir]+1)*i1 + i2-i0 - (i1*(i1+1)/2);
        
        // Set multiplicity
        m = (int(i0)< m_p[dir] ? m_p[dir] - i0 : 0 );
    }

    /*
      Signature
      \a c configuration
      \a m multiplicity
    */
    inline void signature(index_t i0, 
                          index_t i1, 
                          int const & dir, // direction
                          int & c,         // configuration
                          int & m ) const  // multiplicity 
    {
        // Sort indices
        if (i0>i1)
            std::swap(i0,i1);
        
        // Check for right boundary
        if ( int(i1)>m_sz[dir]-m_p[dir] )
        { 
            // mirrior the configuration
            std::swap(i0,i1);
            i0 = m_sz[dir] - i0;
            i1 = m_sz[dir] - i1;
        }
        
        // Set configuration
        c = i1 - i0;
        
        // Set multiplicity
        m = (int(i0)< m_p[dir] ? m_p[dir] - i0 : 0 );
    }

    inline void tri_integral_acc(index_t i0, 
                                  index_t i1, 
                                  index_t i2, 
                                  int const & dir, // direction
                                  gsVector<T> & result ) const
    {
        int c, v, m, q[3];
        q[0]=0;
        q[1]=1;
        q[2]=2;
        bool mirror(false);

        // Sort indices
        if (i0>i1)
        {
            std::swap(i0,i1);
            std::swap(q[0],q[1]);
        }
        if (i1>i2)
        {
            std::swap(i1,i2);
            std::swap(q[1],q[2]);
        }
        if (i0>i1)
        {
            std::swap(i0,i1);
            std::swap(q[0],q[1]);
        }
        
        // Check for right boundary
        if ( int(i2)>m_sz[dir]-m_p[dir] )
        {
            // mirrior the configuration
            std::swap(i0,i2);
            std::swap(q[0],q[2]);
            i0 = m_sz[dir] - i0;
            i1 = m_sz[dir] - i1 - i0;// i1 gets the difference
            i2 = m_sz[dir] - i2;

            mirror = true;
        }
        else
        {
            i1 = i1 - i0;// i1 gets the difference
        }
        
        // Set configuration
        c = (m_p[dir]+1)*i1 + i2-i0 - (i1*(i1+1)/2);
        
        // Set multiplicity
        m = (int(i0)< m_p[dir] ? m_p[dir] - i0 : 0 );

        // ----- Get gradient matrix ----- 
        T corr;

        int qi[3];
        qi[q[0]]=0;
        qi[q[1]]=1;
        qi[q[2]]=2;

        T corr0;
        int v0;
        for ( int r = 0; r!=m_d; ++r)
        {
            corr0 = m_h[dir];
            v0 = 0;
            if ( r == dir )
            {
                v0|= 1 << (2-qi[1]);
                corr0 /= ( mirror ? -m_h[dir] : m_h[dir] );
            }

            for ( int s = 0; s!=m_d; ++s)
            {
                if ( s == dir )
                {
                    v = v0 | (1<< (2-qi[2]) );
                    corr = ( mirror ? -corr0 : corr0 );
                    corr /= m_h[dir];
                    result[ r*m_d + s] *= corr * m_lookup[7*dir+v](c,m) ;
                }
                else
                {
                    result[ r*m_d + s] *= corr0 * m_lookup[7*dir+v0](c,m) ;
                }
            }
        }
    }

    inline void tri_integral(index_t i0, 
                              index_t i1, 
                              index_t i2, 
                              int const & dir, // direction
                              typename gsMatrix<T>::ColXpr result ) const
    {
        //GISMO_ASSERT( result.size() == m_d*m_d; "size not correct.");
        int c, v, m, q[3];
        q[0]=0;
        q[1]=1;
        q[2]=2;
        bool mirror(false);

        // Sort indices
        if (i0>i1)
        {
            std::swap(i0,i1);
            std::swap(q[0],q[1]);
        }
        if (i1>i2)
        {
            std::swap(i1,i2);
            std::swap(q[1],q[2]);
        }
        if (i0>i1)
        {
            std::swap(i0,i1);
            std::swap(q[0],q[1]);
        }
        
        // Check for right boundary
        if ( int(i2)>m_sz[dir]-m_p[dir] )
        {
            // mirrior the configuration
            std::swap(i0,i2);
            std::swap(q[0],q[2]);
            i0 = m_sz[dir] - i0;
            i1 = m_sz[dir] - i1 - i0;// i1 gets the difference
            i2 = m_sz[dir] - i2;

            mirror = true;
        }
        else
        {
            i1 = i1 - i0;// i1 gets the difference
        }
        
        // Set configuration
        c = (m_p[dir]+1)*i1 + i2-i0 - (i1*(i1+1)/2);
        
        // Set multiplicity
        m = (int(i0)< m_p[dir] ? m_p[dir] - i0 : 0 );

        // ----- Get gradient ----- 
        T corr;

        int qi[3];
        qi[q[0]]=0;
        qi[q[1]]=1;
        qi[q[2]]=2;

        T corr0;
        int v0;
        for ( int r = 0; r!=m_d; ++r)
        {
            corr0 = m_h[dir];
            v0 = 0;
            if ( r == dir )
            {
                v0|= 1 << (2- qi[1]);
                corr0 /= ( mirror ? -m_h[dir] : m_h[dir] );
            }

            for ( int s = 0; s!=m_d; ++s)
            {
                if ( s == dir )
                {
                    //v = v0;
                    v= v0 | (1<< (2-qi[2]) );
                    corr = ( mirror ? -corr0 : corr0 );
                    corr /= m_h[dir];
                    result( r*m_d + s,0) = corr * m_lookup[7*dir+v](c,m) ;
                }
                else
                {
                    result( r*m_d + s,0) = corr0 * m_lookup[7*dir+v0](c,m) ;
                }
            }
        }
    }

    inline void normalForm(index_t i0, 
                           index_t i1, 
                           index_t i2, 
                           int const & dir, // direction
                           int & c,
                           int & m,
                           int qi[3]
        ) const
    {
        int q[3];
        q[0]=0;
        q[1]=1;
        q[2]=2;

        // Sort indices
        if (i0>i1)
        {
            std::swap(i0,i1);
            std::swap(q[0],q[1]);
        }
        if (i1>i2)
        {
            std::swap(i1,i2);
            std::swap(q[1],q[2]);
        }
        if (i0>i1)
        {
            std::swap(i0,i1);
            std::swap(q[0],q[1]);
        }
        
        // Check for right boundary
        if ( int(i2)>m_sz[dir]-m_p[dir] )
        {
            // mirrior the configuration
            std::swap(i0,i2);
            std::swap(q[0],q[2]);
            i0 = m_sz[dir] - i0;
            i1 = m_sz[dir] - i1 - i0;// i1 gets the difference
            i2 = m_sz[dir] - i2;
        }
        else
        {
            i1 = i1 - i0;// i1 gets the difference
        }
        
        // Set configuration
        c = (m_p[dir]+1)*i1 + i2-i0 - (i1*(i1+1)/2);
        
        // Set multiplicity
        m = (int(i0)< m_p[dir] ? m_p[dir] - i0 : 0 );

        // Compute permutation
        qi[q[0]]=0;
        qi[q[1]]=1;
        qi[q[2]]=2;
    }

    
    /// Returns multi-index range [K0,K1] overlapping II intersect. JJ
    template<int D> inline
    bool overlapRange(gsVector<index_t,D> const & II, 
                      gsVector<index_t,D> const & JJ, 
                      gsVector<index_t,D> & K0,
                      gsVector<index_t,D> & K1 ) const
    {  
        index_t loffset, roffset;
        
        for (index_t i = 0; i!=II.size(); ++i)
        {
            if ( std::abs( static_cast<int>(II[i]-JJ[i]) ) > m_p[i] )
                return false;
            
            if ( II[i] < JJ[i] )
            {
                loffset = JJ[i]-m_p[i];
                roffset = II[i]+m_p[i];
            }
            else
            {
                loffset = II[i]-m_p[i];
                roffset = JJ[i]+m_p[i];
            }
        
            if (loffset>roffset)
                return false;
        
            K0[i] = ( loffset > 0 ? loffset : 0 );
            K1[i] = ( roffset <= (m_sz[i]) ? roffset : m_sz[i] );	
        }
    
        return true;
    }

    /// Returns multi-index range [K0,K1] overlapping II
    template<int D> inline
    void overlapRange(gsVector<index_t,D> const & II, 
                      gsVector<index_t,D> & K0,
                      gsVector<index_t,D> & K1 ) const
    {  
        for (index_t i = 0; i!=II.size(); ++i)
        {
            const index_t loffset= II[i]-m_p[i], 
                          roffset= II[i]+m_p[i];
            K0[i] = ( loffset >  0         ? loffset : 0       );
            K1[i] = ( roffset <= (m_sz[i]) ? roffset : m_sz[i] );	
        }
    }

    template<int D> inline
    int tindex(gsVector<index_t,D> const & v)
    {
        index_t ind = v(m_d-1);
        for ( index_t i=m_d-2; i>=0; --i )
            ind = ind * (m_sz(i)+1) + v(i) ;
        return ind;
    }

    template<int D> inline
    int isRegular(const gsVector<unsigned,D> & I,
                  const gsVector<unsigned,D> & J,
                  const gsVector<unsigned,D> & K)
    {
        for ( index_t r=0; r != m_d; ++r )
        {
            if ( I[r] < m_p[r] ||
                 J[r] < m_p[r] ||
                 K[r] < m_p[r] ||
                 I[r] > m_sz[r] - m_p[r] ||
                 J[r] > m_sz[r] - m_p[r] ||
                 K[r] > m_sz[r] - m_p[r] )
                return false;
        }
        return true;
    }

    void initialize(const gsBasis<T>& m_basis);

    //bool accept(const gsBasis<T>& m_basis);

public:
  
    gsGeometry<T> * rhsApprox() { return m_rhsApprox; }
    gsGeometry<T> * geoApprox() { return m_geometryFactor; }

// Data members
private:

    /// Pointers to Approximated quantities
    mutable gsGeometry<T> * m_geometryFactor;
    mutable gsGeometry<T> * m_rhsApprox;   
  
    /// Coordinate-wise degrees
    mutable gsVector<int> m_p;
  
    /// Dimension
    index_t m_d;

    /// Coordinate-wise knot-spacing
    mutable gsVector<T> m_h;
  
    /// Coordinate-wise max-indices
    /// (ie. size-1)
    mutable gsVector<index_t> m_sz;
  
    /// Lookup tables
    gsMatrix<T> m_lookup[21]; //    7*m_d 
    gsMatrix<T> m_rhs_lookup[3]; //   m_d
 
    bool m_initialized;
  
}; // class gsInterpolationAssembler


//////////////////////////////////////////////////
//////////////////////////////////////////////////

} // namespace gismo

#ifndef GISMO_HEADERS_ONLY
#include GISMO_HPP_HEADER(gsInterpolationAssembler.hpp)
#endif
